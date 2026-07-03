import contextlib
import io
import json
import os
import sqlite3
import tarfile
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from sync_backend import (
    ModelSettings,
    create_backup,
    load_model_settings,
    main,
    resolve_paths,
    restore_backup,
    sync_history,
    update_rollout_file,
)


class SyncBackendTests(unittest.TestCase):
    def test_main_help_command_prints_global_help(self):
        output = io.StringIO()

        with contextlib.redirect_stdout(output):
            exit_code = main(["help"])

        self.assertEqual(exit_code, 0)
        help_text = output.getvalue()
        self.assertIn("usage: codex-history-sync", help_text)
        self.assertIn("status", help_text)
        self.assertIn("sync", help_text)
        self.assertIn("restore", help_text)

    def test_resolve_paths_expands_environment_variables(self):
        with TemporaryDirectory() as tmp:
            previous = os.environ.get("CODEX_TEST_HOME")
            os.environ["CODEX_TEST_HOME"] = tmp
            try:
                paths = resolve_paths("$CODEX_TEST_HOME/.codex")
            finally:
                if previous is None:
                    os.environ.pop("CODEX_TEST_HOME", None)
                else:
                    os.environ["CODEX_TEST_HOME"] = previous

            self.assertEqual(paths.home, Path(tmp) / ".codex")
            self.assertEqual(paths.config, Path(tmp) / ".codex" / "config.toml")

    def test_load_model_settings_from_config(self):
        with TemporaryDirectory() as tmp:
            config = Path(tmp) / "config.toml"
            config.write_text('model_provider = "anthropic"\nmodel = "claude-sonnet-4-5"\n')

            settings = load_model_settings(config)

            self.assertEqual(settings.provider, "anthropic")
            self.assertEqual(settings.model, "claude-sonnet-4-5")

    def test_sync_updates_database_rollout_and_index(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.config.write_text('model_provider = "openai"\nmodel = "gpt-5.1-codex"\n')
            paths.sessions_dir.mkdir()
            day_dir = paths.sessions_dir / "2026" / "05" / "13"
            day_dir.mkdir(parents=True)
            rollout = day_dir / "rollout-test.jsonl"
            rollout.write_text(
                json.dumps(
                    {
                        "type": "session_meta",
                        "id": "session-1",
                        "payload": {
                            "model_provider": "old",
                            "model": "old-model",
                        },
                    }
                )
                + "\n"
                + json.dumps({"type": "event", "message": "keep me"})
                + "\n"
            )
            paths.session_index.write_text(
                json.dumps(
                    {
                        "id": "session-1",
                        "model_provider": "old",
                        "model": "old-model",
                    }
                )
                + "\n"
            )

            connection = sqlite3.connect(paths.state_db)
            connection.execute(
                "CREATE TABLE threads (id TEXT PRIMARY KEY, model_provider TEXT, model TEXT)"
            )
            connection.execute("INSERT INTO threads VALUES ('session-1', 'old', 'old-model')")
            connection.commit()
            connection.close()

            stats = sync_history(paths, load_model_settings(paths.config))

            self.assertEqual(stats.db_threads_updated, 1)
            self.assertEqual(stats.rollout_files_updated, 1)
            self.assertEqual(stats.index_rows_updated, 1)
            self.assertIsNotNone(stats.backup_path)
            self.assertTrue(stats.backup_path.exists())

            connection = sqlite3.connect(paths.state_db)
            row = connection.execute(
                "SELECT model_provider, model FROM threads WHERE id = 'session-1'"
            ).fetchone()
            connection.close()
            self.assertEqual(row, ("openai", "gpt-5.1-codex"))

            rollout_meta = json.loads(rollout.read_text().splitlines()[0])
            self.assertEqual(rollout_meta["payload"]["model_provider"], "openai")
            self.assertEqual(rollout_meta["payload"]["model"], "gpt-5.1-codex")

            index_meta = json.loads(paths.session_index.read_text().splitlines()[0])
            self.assertEqual(index_meta["model_provider"], "openai")
            self.assertEqual(index_meta["model"], "gpt-5.1-codex")

    def test_dry_run_does_not_modify_files(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.sessions_dir.mkdir()
            rollout = paths.sessions_dir / "rollout-test.jsonl"
            original = json.dumps({"type": "session_meta", "model_provider": "old", "model": "old"})
            rollout.write_text(original + "\n")

            stats = sync_history(
                paths,
                ModelSettings(provider="openai", model="gpt-5"),
                dry_run=True,
            )

            self.assertEqual(stats.rollout_files_updated, 1)
            self.assertEqual(rollout.read_text(), original + "\n")
            self.assertIsNone(stats.backup_path)

    def test_sync_creates_missing_session_index_from_database(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.config.write_text('model_provider = "openai"\nmodel = "gpt-5.5"\n')

            connection = sqlite3.connect(paths.state_db)
            connection.execute(
                "CREATE TABLE threads (id TEXT PRIMARY KEY, title TEXT, updated_at INTEGER, archived INTEGER, model_provider TEXT, model TEXT)"
            )
            connection.execute(
                "INSERT INTO threads VALUES ('session-1', 'Test Session', 1760000000, 0, 'old', 'old-model')"
            )
            connection.commit()
            connection.close()

            self.assertFalse(paths.session_index.exists())

            stats = sync_history(paths, load_model_settings(paths.config))

            self.assertEqual(stats.index_rows_updated, 1)
            self.assertTrue(paths.session_index.exists())

            index_meta = json.loads(paths.session_index.read_text().splitlines()[0])
            self.assertEqual(index_meta["id"], "session-1")
            self.assertEqual(index_meta["thread_name"], "Test Session")

    def test_sync_preserves_git_metadata_in_generated_session_index(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.config.write_text('model_provider = "openai"\nmodel = "gpt-5.5"\n')

            connection = sqlite3.connect(paths.state_db)
            connection.execute(
                "CREATE TABLE threads (id TEXT PRIMARY KEY, title TEXT, updated_at INTEGER, archived INTEGER, model_provider TEXT, model TEXT, cwd TEXT, git_branch TEXT, git_sha TEXT, git_origin_url TEXT, rollout_path TEXT)"
            )
            connection.execute(
                "INSERT INTO threads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    "session-1",
                    "Git Session",
                    1760000000,
                    0,
                    "old",
                    "old-model",
                    "/repo/codex-history-sync",
                    "main",
                    "abc123",
                    "git@gitee.com:duke/codex-history-sync.git",
                    "sessions/2026/05/17/rollout-session-1.jsonl",
                ),
            )
            connection.commit()
            connection.close()

            sync_history(paths, load_model_settings(paths.config))

            index_meta = json.loads(paths.session_index.read_text().splitlines()[0])
            self.assertEqual(index_meta["cwd"], "/repo/codex-history-sync")
            self.assertEqual(index_meta["git_branch"], "main")
            self.assertEqual(index_meta["git_sha"], "abc123")
            self.assertEqual(index_meta["git_origin_url"], "git@gitee.com:duke/codex-history-sync.git")
            self.assertEqual(index_meta["rollout_path"], "sessions/2026/05/17/rollout-session-1.jsonl")
            self.assertEqual(
                index_meta["git"],
                {
                    "branch": "main",
                    "commit_hash": "abc123",
                    "repository_url": "git@gitee.com:duke/codex-history-sync.git",
                },
            )

    def test_sync_keeps_existing_session_index_fields_while_adding_git_metadata(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.config.write_text('model_provider = "openai"\nmodel = "gpt-5.5"\n')
            paths.session_index.write_text(
                json.dumps(
                    {
                        "id": "session-1",
                        "thread_name": "Existing Name",
                        "custom_field": "keep-me",
                        "git": {"dirty": True},
                    }
                )
                + "\n"
            )

            connection = sqlite3.connect(paths.state_db)
            connection.execute(
                "CREATE TABLE threads (id TEXT PRIMARY KEY, title TEXT, updated_at INTEGER, archived INTEGER, model_provider TEXT, model TEXT, git_branch TEXT, git_sha TEXT)"
            )
            connection.execute(
                "INSERT INTO threads VALUES ('session-1', 'Database Name', 1760000000, 0, 'old', 'old-model', 'feature/history', 'def456')"
            )
            connection.commit()
            connection.close()

            sync_history(paths, load_model_settings(paths.config))

            index_meta = json.loads(paths.session_index.read_text().splitlines()[0])
            self.assertEqual(index_meta["thread_name"], "Existing Name")
            self.assertEqual(index_meta["custom_field"], "keep-me")
            self.assertEqual(index_meta["git"]["dirty"], True)
            self.assertEqual(index_meta["git"]["branch"], "feature/history")
            self.assertEqual(index_meta["git"]["commit_hash"], "def456")

    def test_sync_does_not_add_empty_git_object_to_session_index(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.config.write_text('model_provider = "openai"\nmodel = "gpt-5.5"\n')

            connection = sqlite3.connect(paths.state_db)
            connection.execute(
                "CREATE TABLE threads (id TEXT PRIMARY KEY, title TEXT, updated_at INTEGER, archived INTEGER, model_provider TEXT, model TEXT, git_branch TEXT, git_sha TEXT, git_origin_url TEXT)"
            )
            connection.execute(
                "INSERT INTO threads VALUES ('session-1', 'No Git Session', 1760000000, 0, 'old', 'old-model', '', '', '')"
            )
            connection.commit()
            connection.close()

            sync_history(paths, load_model_settings(paths.config))

            index_meta = json.loads(paths.session_index.read_text().splitlines()[0])
            self.assertNotIn("git", index_meta)
            self.assertNotIn("git_branch", index_meta)
            self.assertNotIn("git_sha", index_meta)
            self.assertNotIn("git_origin_url", index_meta)

    def test_restore_rejects_path_traversal(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            backup = paths.home / "bad.tar.gz"
            with tarfile.open(backup, "w:gz") as archive:
                outside = paths.home / "outside.txt"
                outside.write_text("bad")
                archive.add(outside, arcname="../outside.txt")

            with self.assertRaisesRegex(Exception, "outside CODEX_HOME"):
                restore_backup(paths, backup)

    def test_restore_backup_restores_regular_files(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            paths.config.write_text('model_provider = "openai"\n')
            backup = create_backup(paths)
            paths.config.write_text('model_provider = "changed"\n')

            restore_backup(paths, backup)

            self.assertEqual(paths.config.read_text(), 'model_provider = "openai"\n')

    def test_backup_uses_posix_archive_names(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            day_dir = paths.sessions_dir / "2026" / "05" / "13"
            day_dir.mkdir(parents=True)
            rollout = day_dir / "rollout-test.jsonl"
            rollout.write_text("{}\n")

            backup = create_backup(paths)

            with tarfile.open(backup, "r:gz") as archive:
                names = archive.getnames()
            self.assertIn("sessions/2026/05/13/rollout-test.jsonl", names)
            self.assertFalse(any("\\" in name for name in names))

    def test_restore_rejects_windows_backslash_path_traversal(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            backup = paths.home / "bad.tar.gz"
            with tarfile.open(backup, "w:gz") as archive:
                outside = paths.home / "outside.txt"
                outside.write_text("bad")
                archive.add(outside, arcname="..\\outside.txt")

            with self.assertRaisesRegex(Exception, "outside CODEX_HOME"):
                restore_backup(paths, backup)

    def test_restore_rejects_windows_drive_like_member(self):
        with TemporaryDirectory() as tmp:
            paths = resolve_paths(tmp)
            paths.home.mkdir(parents=True, exist_ok=True)
            backup = paths.home / "bad.tar.gz"
            with tarfile.open(backup, "w:gz") as archive:
                outside = paths.home / "outside.txt"
                outside.write_text("bad")
                archive.add(outside, arcname="C:/outside.txt")

            with self.assertRaisesRegex(Exception, "outside CODEX_HOME"):
                restore_backup(paths, backup)

    def test_atomic_rollout_write_preserves_lf_newlines(self):
        with TemporaryDirectory() as tmp:
            rollout = Path(tmp) / "rollout-test.jsonl"
            rollout.write_text(
                json.dumps({"type": "session_meta", "model_provider": "old", "model": "old"})
                + "\n",
                newline="\n",
            )

            changed = update_rollout_file(
                rollout,
                ModelSettings(provider="openai", model="gpt-5"),
            )

            self.assertTrue(changed)
            self.assertNotIn(b"\r\n", rollout.read_bytes())
            self.assertIn(b"\n", rollout.read_bytes())


if __name__ == "__main__":
    unittest.main()
