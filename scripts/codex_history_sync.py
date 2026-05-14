#!/usr/bin/env python3
"""Merge Codex conversation history through a user-selected sync folder.

This script intentionally avoids auth/config/device state. It syncs rollout
JSONL files, lightweight JSONL indexes, and optional insert-only thread metadata
from state_5.sqlite.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import shutil
import socket
import sqlite3
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


SAFE_DIRS = ("sessions", "archived_sessions")
INDEX_FILES = ("session_index.jsonl", "history.jsonl")
STATE_TABLES = (
    "threads",
    "thread_spawn_edges",
    "thread_goals",
    "thread_dynamic_tools",
)
STATE_SNAPSHOT = Path("state") / "state_tables.jsonl"
SCRIPT_VERSION = 1


def utc_stamp() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def utc_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


def resolve_path(value: str | Path) -> Path:
    return Path(value).expanduser().resolve(strict=False)


def default_codex_home() -> Path:
    env_home = os.environ.get("CODEX_HOME")
    if env_home:
        return resolve_path(env_home)
    return resolve_path(Path.home() / ".codex")


def qident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def log(args: argparse.Namespace, message: str) -> None:
    if not getattr(args, "quiet", False):
        print(message)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def same_file_content(a: Path, b: Path) -> bool:
    if not a.exists() or not b.exists():
        return False
    if a.stat().st_size != b.stat().st_size:
        return False
    return sha256_file(a) == sha256_file(b)


def is_allowed_rel(rel: Path) -> bool:
    parts = rel.parts
    if not parts:
        return False
    if len(parts) == 1 and parts[0] in INDEX_FILES:
        return True
    return parts[0] in SAFE_DIRS and ".." not in parts


def iter_session_files(root: Path) -> Iterator[Tuple[Path, Path]]:
    for dirname in SAFE_DIRS:
        base = root / dirname
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if path.is_file():
                yield path.relative_to(root), path


class Counters:
    def __init__(self) -> None:
        self.values: Dict[str, int] = OrderedDict()

    def inc(self, key: str, n: int = 1) -> None:
        self.values[key] = self.values.get(key, 0) + n

    def summary(self) -> str:
        if not self.values:
            return "no changes"
        return ", ".join(f"{k}={v}" for k, v in self.values.items())


class BackupManager:
    def __init__(self, codex_home: Path, label: str, dry_run: bool = False) -> None:
        self.codex_home = codex_home
        self.label = label
        self.dry_run = dry_run
        self.root: Optional[Path] = None

    def ensure(self) -> Path:
        if self.root is None:
            self.root = (
                self.codex_home
                / "backups_state"
                / "history-sync"
                / f"{utc_stamp()}-{self.label}"
            )
            if not self.dry_run:
                self.root.mkdir(parents=True, exist_ok=True)
                manifest = {
                    "tool": "codex-history-sync",
                    "created_at": utc_iso(),
                    "label": self.label,
                    "codex_home": str(self.codex_home),
                }
                (self.root / "manifest.json").write_text(
                    json.dumps(manifest, indent=2), encoding="utf-8"
                )
        return self.root

    def backup_file(self, rel: Path) -> None:
        src = self.codex_home / rel
        if not src.exists():
            return
        dst = self.ensure() / "files" / rel
        if self.dry_run:
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists():
            shutil.copy2(src, dst)

    def backup_sqlite(self, db_path: Path) -> None:
        if not db_path.exists():
            return
        dst = self.ensure() / "state_5.sqlite"
        if self.dry_run:
            return
        backup_sqlite_database(db_path, dst)


def backup_sqlite_database(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        source = sqlite3.connect(f"{src.as_uri()}?mode=ro", uri=True)
        target = sqlite3.connect(str(dst))
        try:
            source.backup(target)
        finally:
            target.close()
            source.close()
    except Exception:
        shutil.copy2(src, dst)
        for suffix in ("-wal", "-shm"):
            sidecar = Path(str(src) + suffix)
            if sidecar.exists():
                shutil.copy2(sidecar, Path(str(dst) + suffix))


def copy_with_policy(
    src: Path,
    dst: Path,
    *,
    rel: Path,
    prefer: str,
    backup: Optional[BackupManager],
    counters: Counters,
    dry_run: bool,
) -> None:
    if not src.exists():
        counters.inc("missing_source")
        return
    if dst.exists() and same_file_content(src, dst):
        counters.inc("unchanged")
        return

    action = "copy"
    if dst.exists():
        if prefer == "dest":
            counters.inc("kept_destination")
            return
        if prefer == "newer":
            src_mtime = src.stat().st_mtime_ns
            dst_mtime = dst.stat().st_mtime_ns
            if dst_mtime > src_mtime:
                counters.inc("kept_newer_destination")
                return
            if dst_mtime == src_mtime:
                conflict = dst.with_name(
                    f"{dst.stem}.conflict-{socket.gethostname()}-{utc_stamp()}{dst.suffix}"
                )
                if not dry_run:
                    conflict.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, conflict)
                counters.inc("conflict_copies")
                return
        action = "overwrite"
        if backup is not None:
            backup.backup_file(rel)

    if not dry_run:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    counters.inc(action)


def parse_jsonl_time(obj: object) -> Tuple[int, str]:
    if not isinstance(obj, dict):
        return (0, "")
    for key in ("updated_at", "created_at"):
        value = obj.get(key)
        if isinstance(value, str):
            normalized = value.replace("Z", "+00:00")
            try:
                dt = _dt.datetime.fromisoformat(normalized)
                return (int(dt.timestamp() * 1000), value)
            except ValueError:
                return (0, value)
    value = obj.get("ts")
    if isinstance(value, (int, float)):
        return (int(value), str(value))
    return (0, "")


def jsonl_key(filename: str, obj: object, raw: str) -> Tuple[str, ...]:
    if isinstance(obj, dict):
        if filename == "session_index.jsonl" and obj.get("id") is not None:
            return ("session_index", str(obj["id"]))
        if filename == "history.jsonl":
            if obj.get("session_id") is not None and obj.get("text") is not None:
                return ("history", str(obj["session_id"]), str(obj["text"]))
            if obj.get("ts") is not None and obj.get("text") is not None:
                return ("history", str(obj["ts"]), str(obj["text"]))
        if obj.get("id") is not None:
            return ("id", str(obj["id"]))
    return ("raw", hashlib.sha256(raw.encode("utf-8", errors="replace")).hexdigest())


def read_jsonl(path: Path, filename: str) -> List[Tuple[Tuple[str, ...], int, int, str]]:
    if not path.exists():
        return []
    entries = []
    try:
        lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    except OSError:
        return []
    for order, line in enumerate(lines):
        raw = line.strip()
        if not raw:
            continue
        try:
            obj = json.loads(raw)
        except json.JSONDecodeError:
            obj = None
        ts, _ = parse_jsonl_time(obj)
        entries.append((jsonl_key(filename, obj, raw), ts, order, raw))
    return entries


def merge_jsonl(
    existing: Path,
    incoming: Path,
    output: Path,
    *,
    rel: Path,
    backup: Optional[BackupManager],
    counters: Counters,
    dry_run: bool,
) -> None:
    filename = rel.name
    merged: Dict[Tuple[str, ...], Tuple[int, int, str]] = OrderedDict()
    sequence = 0
    for path in (existing, incoming):
        for key, ts, _order, raw in read_jsonl(path, filename):
            sequence += 1
            old = merged.get(key)
            if old is None or (ts, sequence) >= (old[0], old[1]):
                merged[key] = (ts, sequence, raw)

    new_lines = [
        item[2]
        for item in sorted(
            merged.values(),
            key=lambda value: (value[0] == 0, value[0], value[1]),
        )
    ]
    new_text = ("\n".join(new_lines) + "\n") if new_lines else ""
    old_text = ""
    if output.exists():
        old_text = output.read_text(encoding="utf-8-sig", errors="replace")
    if old_text == new_text:
        counters.inc("jsonl_unchanged")
        return
    if backup is not None:
        backup.backup_file(rel)
    if not dry_run:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(new_text, encoding="utf-8")
    counters.inc("jsonl_merged")


def table_names(con: sqlite3.Connection) -> set[str]:
    return {
        row[0]
        for row in con.execute(
            "select name from sqlite_master where type='table'"
        ).fetchall()
    }


def table_columns(con: sqlite3.Connection, table: str) -> List[str]:
    return [row[1] for row in con.execute(f"pragma table_info({qident(table)})")]


def table_pk(con: sqlite3.Connection, table: str) -> List[str]:
    rows = con.execute(f"pragma table_info({qident(table)})").fetchall()
    keyed = [(row[5], row[1]) for row in rows if row[5]]
    keyed.sort()
    return [name for _pos, name in keyed]


def rollout_relpath(codex_home: Path, raw_path: object) -> Optional[str]:
    if not isinstance(raw_path, str) or not raw_path:
        return None
    try:
        path = resolve_path(raw_path)
        rel = path.relative_to(codex_home)
    except Exception:
        return None
    if is_allowed_rel(rel):
        return rel.as_posix()
    return None


def write_state_snapshot(codex_home: Path, output: Path) -> int:
    db = codex_home / "state_5.sqlite"
    if not db.exists():
        return 0
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    con = sqlite3.connect(f"{db.as_uri()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        names = table_names(con)
        with output.open("w", encoding="utf-8") as f:
            for table in STATE_TABLES:
                if table not in names:
                    continue
                pk = table_pk(con, table)
                if not pk:
                    continue
                for row in con.execute(f"select * from {qident(table)}"):
                    data = dict(row)
                    record = {
                        "schema": SCRIPT_VERSION,
                        "table": table,
                        "pk": pk,
                        "row": data,
                    }
                    if table == "threads":
                        rel = rollout_relpath(codex_home, data.get("rollout_path"))
                        if rel:
                            record["rollout_relpath"] = rel
                    f.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")))
                    f.write("\n")
                    count += 1
    finally:
        con.close()
    return count


def read_state_records(path: Path) -> List[dict]:
    if not path.exists():
        return []
    records = []
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(record, dict) and isinstance(record.get("row"), dict):
            records.append(record)
    return records


def state_key(record: dict) -> Tuple[str, Tuple[object, ...]]:
    table = str(record.get("table", ""))
    row = record.get("row") or {}
    pk = record.get("pk") or []
    if isinstance(pk, list) and pk:
        return (table, tuple(row.get(col) for col in pk))
    raw = json.dumps(row, sort_keys=True, ensure_ascii=False)
    return (table, (hashlib.sha256(raw.encode("utf-8")).hexdigest(),))


def state_updated_at(record: dict) -> int:
    value = (record.get("row") or {}).get("updated_at")
    if isinstance(value, (int, float)):
        return int(value)
    if isinstance(value, str):
        try:
            return int(float(value))
        except ValueError:
            return 0
    return 0


def write_state_records(path: Path, records: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")))
            f.write("\n")


def merge_state_files(existing: Path, incoming: Path, output: Path, dry_run: bool) -> int:
    merged: Dict[Tuple[str, Tuple[object, ...]], dict] = OrderedDict()
    for record in read_state_records(existing):
        merged[state_key(record)] = record
    for record in read_state_records(incoming):
        key = state_key(record)
        old = merged.get(key)
        if old is None or state_updated_at(record) >= state_updated_at(old):
            merged[key] = record
    ordered = sorted(
        merged.values(),
        key=lambda rec: (
            STATE_TABLES.index(rec.get("table"))
            if rec.get("table") in STATE_TABLES
            else 999,
            state_updated_at(rec),
            repr(state_key(rec)),
        ),
    )
    old_text = existing.read_text(encoding="utf-8", errors="replace") if existing.exists() else ""
    new_text = "".join(
        json.dumps(rec, ensure_ascii=False, separators=(",", ":")) + "\n"
        for rec in ordered
    )
    if old_text != new_text and not dry_run:
        write_state_records(output, ordered)
    return len(ordered)


def import_state_metadata(
    codex_home: Path,
    snapshot: Path,
    *,
    state_conflict: str,
    backup: BackupManager,
    counters: Counters,
    dry_run: bool,
) -> None:
    records = read_state_records(snapshot)
    if not records:
        counters.inc("state_no_snapshot")
        return
    db = codex_home / "state_5.sqlite"
    if not db.exists():
        counters.inc("state_db_missing")
        return

    backup.backup_sqlite(db)
    if dry_run:
        counters.inc("state_rows_planned", len(records))
        return

    con = sqlite3.connect(str(db))
    con.row_factory = sqlite3.Row
    con.execute("pragma busy_timeout=5000")
    try:
        names = table_names(con)
        for record in records:
            table = record.get("table")
            row = dict(record.get("row") or {})
            if table not in names or table not in STATE_TABLES:
                counters.inc("state_table_skipped")
                continue
            if table == "threads" and record.get("rollout_relpath"):
                row["rollout_path"] = str(codex_home / Path(record["rollout_relpath"]))
            columns = set(table_columns(con, table))
            pk = [col for col in (record.get("pk") or []) if col in columns and col in row]
            usable_cols = [col for col in row.keys() if col in columns]
            if not pk or not usable_cols:
                counters.inc("state_row_skipped")
                continue
            where = " and ".join(f"{qident(col)} = ?" for col in pk)
            where_values = [row[col] for col in pk]
            existing = con.execute(
                f"select * from {qident(table)} where {where}", where_values
            ).fetchone()
            if existing is None:
                col_sql = ", ".join(qident(col) for col in usable_cols)
                val_sql = ", ".join("?" for _ in usable_cols)
                con.execute(
                    f"insert into {qident(table)} ({col_sql}) values ({val_sql})",
                    [row[col] for col in usable_cols],
                )
                counters.inc("state_inserted")
                continue
            if state_conflict == "use-newer":
                local_updated = existing["updated_at"] if "updated_at" in existing.keys() else 0
                remote_updated = row.get("updated_at") or 0
                if isinstance(remote_updated, str):
                    try:
                        remote_updated = int(float(remote_updated))
                    except ValueError:
                        remote_updated = 0
                if isinstance(local_updated, str):
                    try:
                        local_updated = int(float(local_updated))
                    except ValueError:
                        local_updated = 0
                if int(remote_updated) > int(local_updated):
                    update_cols = [col for col in usable_cols if col not in pk]
                    if update_cols:
                        set_sql = ", ".join(f"{qident(col)} = ?" for col in update_cols)
                        con.execute(
                            f"update {qident(table)} set {set_sql} where {where}",
                            [row[col] for col in update_cols] + where_values,
                        )
                        counters.inc("state_updated")
                        continue
            counters.inc("state_kept_local")
        con.commit()
    finally:
        con.close()


def store_files_root(sync_root: Path) -> Path:
    return sync_root / "files"


def export_history(args: argparse.Namespace) -> Counters:
    codex_home = args.codex_home
    sync_root = args.sync_root
    counters = Counters()
    files_root = store_files_root(sync_root)
    if not args.dry_run:
        files_root.mkdir(parents=True, exist_ok=True)

    for rel, src in iter_session_files(codex_home):
        dst = files_root / rel
        copy_with_policy(
            src,
            dst,
            rel=rel,
            prefer=args.prefer,
            backup=None,
            counters=counters,
            dry_run=args.dry_run,
        )

    for name in INDEX_FILES:
        rel = Path(name)
        local = codex_home / rel
        shared = files_root / rel
        merge_jsonl(
            shared,
            local,
            shared,
            rel=rel,
            backup=None,
            counters=counters,
            dry_run=args.dry_run,
        )

    if args.include_state:
        with tempfile.TemporaryDirectory() as tmpdir:
            snapshot = Path(tmpdir) / "state_tables.jsonl"
            count = write_state_snapshot(codex_home, snapshot)
            if count:
                if args.dry_run:
                    counters.inc("state_rows_planned", count)
                else:
                    shared_state = sync_root / STATE_SNAPSHOT
                    total = merge_state_files(shared_state, snapshot, shared_state, args.dry_run)
                    counters.inc("state_exported", total)
            else:
                counters.inc("state_no_local_rows")

    if not args.dry_run:
        manifest = {
            "tool": "codex-history-sync",
            "version": SCRIPT_VERSION,
            "updated_at": utc_iso(),
            "last_export_device": socket.gethostname(),
            "codex_home": str(codex_home),
            "include_state": bool(args.include_state),
        }
        sync_root.mkdir(parents=True, exist_ok=True)
        (sync_root / "manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )
    return counters


def import_history(args: argparse.Namespace) -> Tuple[Counters, Optional[Path]]:
    codex_home = args.codex_home
    sync_root = args.sync_root
    counters = Counters()
    backup = BackupManager(codex_home, "import", args.dry_run)
    files_root = store_files_root(sync_root)

    for dirname in SAFE_DIRS:
        base = files_root / dirname
        if not base.exists():
            continue
        for src in base.rglob("*"):
            if not src.is_file():
                continue
            rel = src.relative_to(files_root)
            if not is_allowed_rel(rel):
                counters.inc("unsafe_skipped")
                continue
            copy_with_policy(
                src,
                codex_home / rel,
                rel=rel,
                prefer=args.prefer,
                backup=backup,
                counters=counters,
                dry_run=args.dry_run,
            )

    for name in INDEX_FILES:
        rel = Path(name)
        shared = files_root / rel
        local = codex_home / rel
        merge_jsonl(
            local,
            shared,
            local,
            rel=rel,
            backup=backup,
            counters=counters,
            dry_run=args.dry_run,
        )

    if args.include_state:
        state_snapshot = sync_root / STATE_SNAPSHOT
        if args.command == "sync" and args.dry_run and not state_snapshot.exists():
            counters.inc("state_import_skipped_dry_run")
        else:
            import_state_metadata(
                codex_home,
                state_snapshot,
                state_conflict=args.state_conflict,
                backup=backup,
                counters=counters,
                dry_run=args.dry_run,
            )

    return counters, backup.root


def status(args: argparse.Namespace) -> None:
    codex_home = args.codex_home
    sync_root = args.sync_root
    print(f"Codex home: {codex_home}")
    print(f"Sync root:  {sync_root if sync_root else '(not provided)'}")
    print("")
    local_sessions = sum(1 for _ in iter_session_files(codex_home))
    print(f"Local session files: {local_sessions}")
    for name in INDEX_FILES:
        path = codex_home / name
        print(f"Local {name}: {'yes' if path.exists() else 'no'}")
    print(f"Local state_5.sqlite: {'yes' if (codex_home / 'state_5.sqlite').exists() else 'no'}")
    if not sync_root:
        return
    files_root = store_files_root(sync_root)
    remote_sessions = 0
    for dirname in SAFE_DIRS:
        base = files_root / dirname
        if base.exists():
            remote_sessions += sum(1 for p in base.rglob("*") if p.is_file())
    print("")
    print(f"Sync session files: {remote_sessions}")
    for name in INDEX_FILES:
        path = files_root / name
        print(f"Sync {name}: {'yes' if path.exists() else 'no'}")
    state_path = sync_root / STATE_SNAPSHOT
    print(f"Sync state snapshot: {'yes' if state_path.exists() else 'no'}")
    manifest = sync_root / "manifest.json"
    if manifest.exists():
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            print(f"Last export: {data.get('updated_at', '(unknown)')}")
            print(f"Last export device: {data.get('last_export_device', '(unknown)')}")
        except Exception:
            print("Manifest: unreadable")


def restore(args: argparse.Namespace) -> Counters:
    backup_dir = args.backup_dir
    codex_home = args.codex_home
    counters = Counters()
    files = backup_dir / "files"
    if files.exists():
        for src in files.rglob("*"):
            if not src.is_file():
                continue
            rel = src.relative_to(files)
            if not is_allowed_rel(rel):
                counters.inc("unsafe_skipped")
                continue
            dst = codex_home / rel
            if not args.dry_run:
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
            counters.inc("files_restored")
    state_backup = backup_dir / "state_5.sqlite"
    if state_backup.exists() and not args.no_state:
        target = codex_home / "state_5.sqlite"
        if not args.dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            for suffix in ("-wal", "-shm"):
                sidecar = Path(str(target) + suffix)
                if sidecar.exists():
                    sidecar.unlink()
            shutil.copy2(state_backup, target)
        counters.inc("state_restored")
    return counters


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--codex-home",
        type=resolve_path,
        default=default_codex_home(),
        help="Codex home directory; defaults to CODEX_HOME or ~/.codex.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Plan changes without writing.")
    parser.add_argument("--quiet", action="store_true", help="Only print command summaries.")


def add_sync_root(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--sync-root",
        type=resolve_path,
        required=True,
        help="Private shared folder used as the history sync store.",
    )


def add_merge_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--prefer",
        choices=("newer", "source", "dest"),
        default="newer",
        help="File conflict policy. 'source' overwrites, 'dest' keeps existing.",
    )
    parser.add_argument(
        "--include-state",
        action="store_true",
        help="Also sync insert-only state_5.sqlite thread visibility metadata.",
    )
    parser.add_argument(
        "--state-conflict",
        choices=("keep-local", "use-newer"),
        default="keep-local",
        help="How to handle existing SQLite thread metadata rows during import.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Synchronize Codex conversation history through a shared folder."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("status", help="Inspect local and sync-folder history state.")
    add_common(p)
    p.add_argument("--sync-root", type=resolve_path, help="Optional sync store to inspect.")

    p = sub.add_parser("export", help="Merge local Codex history into the sync folder.")
    add_common(p)
    add_sync_root(p)
    add_merge_options(p)

    p = sub.add_parser("import", help="Merge sync-folder history into local Codex history.")
    add_common(p)
    add_sync_root(p)
    add_merge_options(p)

    p = sub.add_parser("sync", help="Export then import for bidirectional merge.")
    add_common(p)
    add_sync_root(p)
    add_merge_options(p)

    p = sub.add_parser("restore", help="Restore from a backup made by this script.")
    add_common(p)
    p.add_argument("--backup-dir", type=resolve_path, required=True)
    p.add_argument("--no-state", action="store_true", help="Do not restore state_5.sqlite.")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.codex_home = resolve_path(args.codex_home)
    if hasattr(args, "sync_root") and args.sync_root:
        args.sync_root = resolve_path(args.sync_root)

    try:
        if args.command == "status":
            status(args)
            return 0
        if args.command == "export":
            counters = export_history(args)
            log(args, f"Export complete: {counters.summary()}")
            return 0
        if args.command == "import":
            counters, backup_path = import_history(args)
            log(args, f"Import complete: {counters.summary()}")
            if backup_path:
                log(args, f"Backup: {backup_path}")
            return 0
        if args.command == "sync":
            export_counters = export_history(args)
            import_counters, backup_path = import_history(args)
            log(args, f"Export complete: {export_counters.summary()}")
            log(args, f"Import complete: {import_counters.summary()}")
            if backup_path:
                log(args, f"Backup: {backup_path}")
            return 0
        if args.command == "restore":
            counters = restore(args)
            log(args, f"Restore complete: {counters.summary()}")
            return 0
    except sqlite3.OperationalError as exc:
        print(f"SQLite error: {exc}", file=sys.stderr)
        print("Close Codex/Codex Desktop and retry if the database is busy.", file=sys.stderr)
        return 2
    except PermissionError as exc:
        print(f"Permission error: {exc}", file=sys.stderr)
        print("Close Codex or the sync client if it is holding the file.", file=sys.stderr)
        return 3
    except OSError as exc:
        print(f"File error: {exc}", file=sys.stderr)
        return 4
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
