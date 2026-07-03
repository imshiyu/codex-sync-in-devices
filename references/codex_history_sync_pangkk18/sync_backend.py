#!/usr/bin/env python3
"""Synchronize Codex history metadata on Windows, Linux, and macOS.

Codex stores the selected model/provider in multiple places:

* ~/.codex/config.toml stores the current defaults.
* ~/.codex/state_5.sqlite stores indexed conversation metadata.
* ~/.codex/sessions/**/rollout-*.jsonl stores per-session metadata.
* ~/.codex/session_index.jsonl stores the global history index.

This command-line tool copies the current config model/provider into the
history metadata so older sessions show the same provider/model consistently.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sqlite3
import sys
import tarfile
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback.
    tomllib = None  # type: ignore[assignment]


DEFAULT_PROVIDER = "openai"
DEFAULT_MODEL = "gpt-5"
STATE_DB_NAME = "state_5.sqlite"
SESSION_INDEX_NAME = "session_index.jsonl"
BACKUP_DIR_NAME = "history-sync-backups"
UTC = timezone.utc


class SyncError(RuntimeError):
    """Raised for recoverable user-facing sync failures."""


@dataclass(frozen=True)
class CodexPaths:
    home: Path
    config: Path
    state_db: Path
    sessions_dir: Path
    session_index: Path
    backup_dir: Path


@dataclass(frozen=True)
class ModelSettings:
    provider: str
    model: str


@dataclass
class SyncStats:
    db_threads_seen: int = 0
    db_threads_updated: int = 0
    rollout_files_seen: int = 0
    rollout_files_updated: int = 0
    index_rows_seen: int = 0
    index_rows_updated: int = 0
    malformed_json_lines: int = 0
    backup_path: Path | None = None

    @property
    def changed(self) -> bool:
        return any(
            (
                self.db_threads_updated,
                self.rollout_files_updated,
                self.index_rows_updated,
            )
        )


def resolve_paths(codex_home: str | None = None) -> CodexPaths:
    home_value = codex_home or os.environ.get("CODEX_HOME", "~/.codex")
    home = Path(os.path.expandvars(home_value)).expanduser()
    return CodexPaths(
        home=home,
        config=home / "config.toml",
        state_db=home / STATE_DB_NAME,
        sessions_dir=home / "sessions",
        session_index=home / SESSION_INDEX_NAME,
        backup_dir=home / BACKUP_DIR_NAME,
    )


def load_model_settings(config_path: Path) -> ModelSettings:
    """Read provider/model from config.toml with conservative defaults."""
    if not config_path.exists():
        return ModelSettings(provider=DEFAULT_PROVIDER, model=DEFAULT_MODEL)

    data = _load_toml(config_path)

    provider = _first_string(
        data,
        (
            ("model_provider",),
            ("modelProvider",),
            ("provider",),
            ("defaults", "model_provider"),
            ("defaults", "provider"),
        ),
    )
    model = _first_string(
        data,
        (
            ("model",),
            ("defaults", "model"),
        ),
    )

    return ModelSettings(
        provider=(provider or DEFAULT_PROVIDER).strip() or DEFAULT_PROVIDER,
        model=(model or DEFAULT_MODEL).strip() or DEFAULT_MODEL,
    )


def _load_toml(config_path: Path) -> dict[str, Any]:
    raw = config_path.read_text(encoding="utf-8")
    if tomllib is not None:
        return tomllib.loads(raw)

    data: dict[str, Any] = {}
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip("\"'")
        if key in {"model", "model_provider", "modelProvider", "provider"}:
            data[key] = value
    return data


def _first_string(data: dict[str, Any], paths: Iterable[tuple[str, ...]]) -> str | None:
    for path in paths:
        current: Any = data
        for key in path:
            if not isinstance(current, dict) or key not in current:
                current = None
                break
            current = current[key]
        if isinstance(current, str):
            return current
    return None


def create_backup(paths: CodexPaths) -> Path:
    """Create a compressed backup of files this tool may modify."""
    paths.backup_dir.mkdir(parents=True, exist_ok=True)
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    backup_path = paths.backup_dir / f"codex-history-{timestamp}.tar.gz"

    with tarfile.open(backup_path, "w:gz") as archive:
        for file_path in _backup_candidates(paths):
            if file_path.exists():
                archive.add(file_path, arcname=_archive_name(paths, file_path))

    return backup_path


def _archive_name(paths: CodexPaths, file_path: Path) -> str:
    return file_path.relative_to(paths.home).as_posix()


def _backup_candidates(paths: CodexPaths) -> Iterable[Path]:
    yield paths.config
    yield paths.state_db
    yield paths.session_index
    if paths.sessions_dir.exists():
        yield from paths.sessions_dir.rglob("rollout-*.jsonl")


def list_backups(paths: CodexPaths) -> list[Path]:
    if not paths.backup_dir.exists():
        return []
    return sorted(paths.backup_dir.glob("codex-history-*.tar.gz"))


def restore_backup(paths: CodexPaths, backup_path: Path) -> None:
    backup_path = Path(os.path.expandvars(str(backup_path))).expanduser()
    if not backup_path.exists():
        raise SyncError(f"Backup not found: {backup_path}")
    if not tarfile.is_tarfile(backup_path):
        raise SyncError(f"Backup is not a tar archive: {backup_path}")

    with tarfile.open(backup_path, "r:gz") as archive:
        for member in archive.getmembers():
            _validate_restore_member(paths, member)
        for member in archive.getmembers():
            target = _restore_target(paths, member.name)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            source = archive.extractfile(member)
            if source is None:
                raise SyncError(f"Unable to read backup member: {member.name}")
            target.parent.mkdir(parents=True, exist_ok=True)
            with source, target.open("wb") as output:
                shutil.copyfileobj(source, output)
            if member.mode:
                try:
                    target.chmod(member.mode)
                except OSError:
                    pass


def _validate_restore_member(paths: CodexPaths, member: tarfile.TarInfo) -> None:
    if not (member.isfile() or member.isdir()):
        raise SyncError(f"Refusing to restore non-regular backup member: {member.name}")
    _restore_target(paths, member.name)


def _restore_target(paths: CodexPaths, member_name: str) -> Path:
    normalized_name = member_name.replace("\\", "/")
    posix_path = PurePosixPath(normalized_name)
    parts = posix_path.parts
    if (
        not normalized_name
        or posix_path.is_absolute()
        or ".." in parts
        or _has_windows_drive_prefix(parts)
    ):
        raise SyncError(f"Refusing to restore path outside CODEX_HOME: {member_name}")

    target = (paths.home.joinpath(*parts)).resolve()
    home = paths.home.resolve()
    if target != home and home not in target.parents:
        raise SyncError(f"Refusing to restore path outside CODEX_HOME: {member_name}")
    return target


def _has_windows_drive_prefix(parts: tuple[str, ...]) -> bool:
    return bool(parts and len(parts[0]) == 2 and parts[0][1] == ":" and parts[0][0].isalpha())


def sync_history(paths: CodexPaths, settings: ModelSettings, *, dry_run: bool = False) -> SyncStats:
    if not paths.home.exists():
        raise SyncError(f"Codex home does not exist: {paths.home}")

    stats = SyncStats()
    if not dry_run:
        stats.backup_path = create_backup(paths)

    sync_state_database(paths, settings, stats, dry_run=dry_run)
    sync_rollout_files(paths, settings, stats, dry_run=dry_run)
    sync_session_index(paths, settings, stats, dry_run=dry_run)
    return stats


def sync_state_database(
    paths: CodexPaths,
    settings: ModelSettings,
    stats: SyncStats,
    *,
    dry_run: bool = False,
) -> None:
    if not paths.state_db.exists():
        return

    connection = sqlite3.connect(paths.state_db, timeout=30.0)
    try:
        connection.execute("PRAGMA busy_timeout = 30000")
        columns = _table_columns(connection, "threads")
        if not {"id", "model_provider", "model"}.issubset(columns):
            return

        rows = connection.execute("SELECT id, model_provider, model FROM threads").fetchall()
        stats.db_threads_seen = len(rows)
        to_update = [
            row_id
            for row_id, provider, model in rows
            if provider != settings.provider or model != settings.model
        ]
        stats.db_threads_updated = len(to_update)

        if to_update and not dry_run:
            connection.execute("BEGIN IMMEDIATE")
            connection.executemany(
                "UPDATE threads SET model_provider = ?, model = ? WHERE id = ?",
                ((settings.provider, settings.model, row_id) for row_id in to_update),
            )
            connection.commit()
    finally:
        connection.close()


def _table_columns(connection: sqlite3.Connection, table_name: str) -> set[str]:
    rows = connection.execute(f"PRAGMA table_info({table_name})").fetchall()
    return {str(row[1]) for row in rows}


def sync_rollout_files(
    paths: CodexPaths,
    settings: ModelSettings,
    stats: SyncStats,
    *,
    dry_run: bool = False,
) -> None:
    if not paths.sessions_dir.exists():
        return

    for file_path in sorted(paths.sessions_dir.rglob("rollout-*.jsonl")):
        stats.rollout_files_seen += 1
        changed = update_rollout_file(file_path, settings, dry_run=dry_run)
        if changed:
            stats.rollout_files_updated += 1


def update_rollout_file(file_path: Path, settings: ModelSettings, *, dry_run: bool = False) -> bool:
    lines = file_path.read_text(encoding="utf-8").splitlines(keepends=True)
    if not lines:
        return False

    try:
        first_record = json.loads(lines[0])
    except json.JSONDecodeError:
        return False

    if not isinstance(first_record, dict) or first_record.get("type") != "session_meta":
        return False

    payload = first_record.get("payload")
    target = payload if isinstance(payload, dict) else first_record
    changed = _apply_model_fields(target, settings, add_missing=True)
    if changed and not dry_run:
        lines[0] = json.dumps(first_record, ensure_ascii=False, separators=(",", ":")) + "\n"
        _atomic_write_text(file_path, "".join(lines))
    return changed


def sync_session_index(
    paths: CodexPaths,
    settings: ModelSettings,
    stats: SyncStats,
    *,
    dry_run: bool = False,
) -> None:
    if not paths.session_index.exists() and not paths.state_db.exists():
        return

    existing_lines = paths.session_index.read_text(encoding="utf-8").splitlines() if paths.session_index.exists() else []
    existing_entries: dict[str, dict[str, Any]] = {}
    existing_order: list[str] = []

    for line in existing_lines:
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            stats.malformed_json_lines += 1
            continue
        if not isinstance(record, dict):
            continue
        thread_id = str(record.get("id") or "").strip()
        if not thread_id:
            continue
        stats.index_rows_seen += 1
        existing_entries[thread_id] = record
        existing_order.append(thread_id)

    db_entries = _read_index_entries_from_database(paths, existing_entries)
    if db_entries is None:
        output = []
        for thread_id in existing_order:
            record = dict(existing_entries[thread_id])
            _apply_model_fields(record, settings, add_missing=False)
            output.append(record)
    else:
        db_ids = {str(entry["id"]) for entry in db_entries}
        index_only_entries = [
            existing_entries[thread_id]
            for thread_id in existing_order
            if thread_id not in db_ids
        ]
        output = db_entries + index_only_entries
        for entry in output:
            _apply_model_fields(entry, settings, add_missing=False)
        output.sort(key=lambda item: (_parse_index_timestamp(str(item.get("updated_at") or "")), str(item.get("id") or "")))

    current_text = "\n".join(existing_lines)
    desired_lines = [json.dumps(entry, ensure_ascii=False, separators=(",", ":")) for entry in output]
    desired_text = "\n".join(desired_lines)
    if desired_text:
        desired_text += "\n"
    if current_text and not current_text.endswith("\n"):
        current_text += "\n"

    if desired_text != current_text:
        stats.index_rows_updated = len(output)
        if not dry_run:
            _atomic_write_text(paths.session_index, desired_text)


def _read_index_entries_from_database(
    paths: CodexPaths,
    existing_entries: dict[str, dict[str, Any]],
) -> list[dict[str, Any]] | None:
    if not paths.state_db.exists():
        return None

    connection = sqlite3.connect(paths.state_db, timeout=30.0)
    try:
        connection.row_factory = sqlite3.Row
        columns = _table_columns(connection, "threads")
        if "id" not in columns:
            return None

        selected = ["id"]
        if "title" in columns:
            selected.append("title")
        if "updated_at" in columns:
            selected.append("updated_at")
        for column in ("cwd", "git_branch", "git_sha", "git_origin_url", "rollout_path"):
            if column in columns:
                selected.append(column)
        where_sql = "WHERE archived = 0" if "archived" in columns else ""
        rows = connection.execute(
            f"SELECT {', '.join(selected)} FROM threads {where_sql} ORDER BY id ASC"
        ).fetchall()
    finally:
        connection.close()

    entries: list[dict[str, Any]] = []
    for row in rows:
        thread_id = str(row["id"])
        existing = dict(existing_entries.get(thread_id) or {})
        title = str(row["title"]) if "title" in row.keys() and row["title"] else thread_id
        updated_at = (
            _iso_utc_from_unix(row["updated_at"])
            if "updated_at" in row.keys() and row["updated_at"]
            else str(existing.get("updated_at") or "")
        )
        existing["id"] = thread_id
        existing["thread_name"] = str(existing.get("thread_name") or title)
        existing["updated_at"] = updated_at
        _apply_thread_metadata(existing, row)
        entries.append(existing)
    return entries


def _apply_thread_metadata(entry: dict[str, Any], row: sqlite3.Row) -> None:
    row_keys = set(row.keys())
    for key in ("cwd", "git_branch", "git_sha", "git_origin_url", "rollout_path"):
        if key in row_keys and row[key]:
            entry[key] = str(row[key])

    git_metadata = {
        "branch": entry.get("git_branch"),
        "commit_hash": entry.get("git_sha"),
        "repository_url": entry.get("git_origin_url"),
    }
    git_metadata = {
        key: str(value)
        for key, value in git_metadata.items()
        if value
    }
    if git_metadata:
        existing_git = entry.get("git")
        git = dict(existing_git) if isinstance(existing_git, dict) else {}
        git.update(git_metadata)
        entry["git"] = git


def _iso_utc_from_unix(value: Any) -> str:
    timestamp = int(value)
    if timestamp > 10_000_000_000:
        timestamp = timestamp // 1000
    return datetime.fromtimestamp(timestamp, tz=UTC).isoformat().replace("+00:00", "Z")


def _parse_index_timestamp(value: str) -> datetime:
    if not value:
        return datetime.fromtimestamp(0, tz=UTC)
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return datetime.fromtimestamp(0, tz=UTC)
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=UTC)
    return parsed.astimezone(UTC)


def _apply_model_fields(
    record: dict[str, Any],
    settings: ModelSettings,
    *,
    add_missing: bool,
) -> bool:
    changed = False
    for provider_key in ("model_provider", "modelProvider", "provider"):
        if provider_key in record and record.get(provider_key) != settings.provider:
            record[provider_key] = settings.provider
            changed = True
    for model_key in ("model", "model_name", "modelName"):
        if model_key in record and record.get(model_key) != settings.model:
            record[model_key] = settings.model
            changed = True

    if add_missing and not any(key in record for key in ("model_provider", "modelProvider", "provider")):
        record["model_provider"] = settings.provider
        changed = True
    if add_missing and not any(key in record for key in ("model", "model_name", "modelName")):
        record["model"] = settings.model
        changed = True
    return changed


def _atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        if path.exists():
            shutil.copystat(path, temp_path, follow_symlinks=False)
        else:
            temp_path.chmod(0o644)
        temp_path.replace(path)
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="codex-history-sync",
        description="Synchronize Codex history model/provider metadata on Windows, Linux, and macOS.",
    )
    parser.add_argument(
        "--codex-home",
        help="Codex data directory. Defaults to CODEX_HOME or ~/.codex.",
    )

    subparsers = parser.add_subparsers(dest="command")

    status_parser = subparsers.add_parser("status", help="Show detected paths and current settings.")
    status_parser.set_defaults(func=command_status)

    sync_parser = subparsers.add_parser("sync", help="Synchronize history metadata.")
    sync_parser.add_argument("--dry-run", action="store_true", help="Show what would change without writing.")
    sync_parser.set_defaults(func=command_sync)

    backup_parser = subparsers.add_parser("backup", help="Create a backup without modifying history.")
    backup_parser.set_defaults(func=command_backup)

    list_parser = subparsers.add_parser("list-backups", help="List available backups.")
    list_parser.set_defaults(func=command_list_backups)

    restore_parser = subparsers.add_parser("restore", help="Restore a backup archive.")
    restore_parser.add_argument("backup", help="Path to a backup created by this tool.")
    restore_parser.set_defaults(func=command_restore)

    help_parser = subparsers.add_parser("help", help="Show this help message and exit.")
    help_parser.set_defaults(func=command_help, parser=parser)

    parser.set_defaults(func=command_sync)
    return parser


def command_help(args: argparse.Namespace) -> int:
    args.parser.print_help()
    return 0


def command_status(args: argparse.Namespace) -> int:
    paths = resolve_paths(args.codex_home)
    settings = load_model_settings(paths.config)
    print(f"Codex home:      {paths.home}")
    print(f"Config:          {paths.config} ({_exists(paths.config)})")
    print(f"State database:  {paths.state_db} ({_exists(paths.state_db)})")
    print(f"Sessions dir:    {paths.sessions_dir} ({_exists(paths.sessions_dir)})")
    print(f"Session index:   {paths.session_index} ({_exists(paths.session_index)})")
    print(f"Model provider:  {settings.provider}")
    print(f"Model:           {settings.model}")
    return 0


def command_sync(args: argparse.Namespace) -> int:
    paths = resolve_paths(args.codex_home)
    settings = load_model_settings(paths.config)
    dry_run = bool(getattr(args, "dry_run", False))
    stats = sync_history(paths, settings, dry_run=dry_run)

    mode = "Dry run" if dry_run else "Synced"
    print(f"{mode}: provider={settings.provider}, model={settings.model}")
    print(f"Database threads: {stats.db_threads_updated}/{stats.db_threads_seen} updated")
    print(f"Rollout files:    {stats.rollout_files_updated}/{stats.rollout_files_seen} updated")
    print(f"Index rows:       {stats.index_rows_updated}/{stats.index_rows_seen} updated")
    if stats.malformed_json_lines:
        print(f"Malformed index lines skipped: {stats.malformed_json_lines}")
    if stats.backup_path:
        print(f"Backup: {stats.backup_path}")
    if dry_run and not stats.changed:
        print("No changes needed.")
    return 0


def command_backup(args: argparse.Namespace) -> int:
    paths = resolve_paths(args.codex_home)
    if not paths.home.exists():
        raise SyncError(f"Codex home does not exist: {paths.home}")
    backup_path = create_backup(paths)
    print(f"Backup: {backup_path}")
    return 0


def command_list_backups(args: argparse.Namespace) -> int:
    paths = resolve_paths(args.codex_home)
    backups = list_backups(paths)
    if not backups:
        print("No backups found.")
        return 0
    for backup in backups:
        print(backup)
    return 0


def command_restore(args: argparse.Namespace) -> int:
    paths = resolve_paths(args.codex_home)
    restore_backup(paths, Path(args.backup))
    print(f"Restored: {args.backup}")
    return 0


def _exists(path: Path) -> str:
    return "found" if path.exists() else "missing"


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except SyncError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    except sqlite3.DatabaseError as error:
        print(f"SQLite error: {error}", file=sys.stderr)
        return 3
    except OSError as error:
        print(f"I/O error: {error}", file=sys.stderr)
        return 4


if __name__ == "__main__":
    raise SystemExit(main())
