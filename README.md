# Codex History Sync

Codex skill for synchronizing conversation history across devices through a private shared folder such as OneDrive, Dropbox, Syncthing, iCloud Drive, a private Git repo, or an external drive.

The skill syncs Codex conversation history and optional thread visibility metadata. It intentionally excludes login, account, device, provider, and configuration state.

## Install

Copy this folder to your Codex skills directory:

```powershell
Copy-Item -Recurse -Force . "$env:USERPROFILE\.codex\skills\codex-history-sync"
```

On macOS/Linux:

```bash
cp -R . ~/.codex/skills/codex-history-sync
```

Restart Codex or open a new session so the skill list is refreshed.

## Quick Usage

Pick a private shared folder, for example `D:\Sync\codex-history`.

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" status --sync-root "D:\Sync\codex-history"
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" export --sync-root "D:\Sync\codex-history"
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" import --sync-root "D:\Sync\codex-history"
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history"
```

If imported sessions do not appear in Codex history or `/resume`, include thread metadata:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history" --include-state
```

Use `--dry-run` before a first import:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history" --dry-run
```

## Safety Boundary

Synced by default:

- `sessions/**`
- `archived_sessions/**`
- `session_index.jsonl`
- `history.jsonl`

Optional with `--include-state`:

- insert-only selected thread metadata from `state_5.sqlite`

Never synced:

- `auth.json`
- `config.toml`
- `installation_id`
- `cap_sid`
- `.sandbox*`
- `cache`
- `log`
- `logs_*.sqlite`
- `tmp`
- `skills`
- `memories`

Imports create backups under:

```text
~/.codex/backups_state/history-sync/
```

## Credits

Inspired by the safety boundary and sync goals of [Dailin521/codex-provider-sync](https://github.com/Dailin521/codex-provider-sync).

