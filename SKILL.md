---
name: codex-history-sync
description: Use when the user wants to synchronize, back up, import, export, restore, or diagnose Codex conversation history across devices or cloud folders. Handles ~/.codex/sessions, archived_sessions, session_index.jsonl, history.jsonl, and optional insert-only state_5.sqlite thread metadata while excluding auth.json, config.toml, installation_id, provider credentials, and other device/login state.
---

# Codex History Sync

## Overview

Synchronize Codex conversation history through a user-chosen folder such as OneDrive, Dropbox, iCloud Drive, Syncthing, a private Git repo, or an external drive. The bundled script performs merge-style export/import and creates local backups before import-side writes.

This skill is inspired by the safety boundary used by `Dailin521/codex-provider-sync`: treat rollout files and thread metadata as syncable history, but do not touch auth, account, provider, or encrypted-content ownership state.

## Safety Rules

- Never sync or edit `auth.json`, `config.toml`, `installation_id`, `cap_sid`, `.sandbox*`, `cache`, `log`, `tmp`, `skills`, or `memories`.
- Default sync includes rollout JSONL files and lightweight JSONL indexes only.
- Use `--include-state` only when the user wants imported conversations to appear in Codex Desktop or `/resume` and accepts SQLite metadata writes.
- State metadata import is insert-only by default. It adds missing thread rows and does not overwrite existing local thread rows unless `--state-conflict use-newer` is explicitly selected.
- Do not use this skill to move histories across different accounts/providers when the user expects encrypted content to become decryptable. Old rows containing provider/account-bound `encrypted_content` can remain visible but may not be resumable.
- If Codex or Codex Desktop is actively writing history, prefer closing it before `import`, `sync --include-state`, or `restore`.

## Workflow

1. Ask for or infer the shared sync folder. It should be a private directory controlled by the user.
2. Run `status` first to inspect the local Codex home and the sync folder.
3. On the first device, run `export`.
4. On other devices, run `import`, or run `sync` for bidirectional merge.
5. Add `--include-state` only if history files import correctly but conversations are not visible in Codex UI/history lists.
6. If an import behaves badly, use the printed backup path with `restore`.

## Commands

Use the bundled script:

```powershell
python .\scripts\codex_history_sync.py status --sync-root "D:\Sync\codex-history"
python .\scripts\codex_history_sync.py export --sync-root "D:\Sync\codex-history"
python .\scripts\codex_history_sync.py import --sync-root "D:\Sync\codex-history"
python .\scripts\codex_history_sync.py sync --sync-root "D:\Sync\codex-history"
```

For Codex Desktop visibility metadata:

```powershell
python .\scripts\codex_history_sync.py sync --sync-root "D:\Sync\codex-history" --include-state
```

For dry runs:

```powershell
python .\scripts\codex_history_sync.py import --sync-root "D:\Sync\codex-history" --dry-run
```

For restore:

```powershell
python .\scripts\codex_history_sync.py restore --backup-dir "$env:USERPROFILE\.codex\backups_state\history-sync\<timestamp>"
```

On macOS/Linux, use the same commands with `python3` and POSIX paths.

## Script Behavior

- `export`: merges local history into the sync folder without deleting remote files.
- `import`: merges the sync folder into local Codex history and backs up overwritten index/session files.
- `sync`: runs export then import for a bidirectional union.
- `status`: reports local and sync-folder counts.
- `restore`: restores files from a backup created by this script.

The script only depends on Python standard library modules.

## References

Read `references/safety.md` before using `--include-state`, resolving conflicts, or restoring backups.
