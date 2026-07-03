# Codex History Sync

English | [简体中文](README.md)

A command-line Codex history synchronization tool for Windows, Linux, Ubuntu, and macOS. It has no UI and is designed for local developer machines, servers, and scriptable environments.

This tool synchronizes the current Codex model metadata into historical session metadata. It helps fix stale or inconsistent `model_provider` / `model` values shown in Codex history files, indexes, and session records.

Codex history metadata is usually stored in:

- `~/.codex/config.toml`
- `~/.codex/state_5.sqlite`
- `~/.codex/sessions/**/rollout-*.jsonl`
- `~/.codex/session_index.jsonl`

Before any real write, the tool creates a `tar.gz` backup automatically.

## Features

- Supports Windows, Linux, Ubuntu, and macOS.
- No third-party dependencies; Python is enough.
- Reads `CODEX_HOME` automatically, falling back to `~/.codex`.
- Synchronizes the SQLite history database, session JSONL files, and the global session index.
- Preserves historical context such as working directory, Git branch, commit hash, remote URL, and rollout path when rebuilding the session index.
- Supports `--dry-run`.
- Creates automatic backups before modifying files.
- Supports listing and restoring backups.

## Requirements

- Python 3.9+
- No third-party Python packages

## Get the Code

From GitHub:

```bash
git clone git@github.com:pangkk18/codex-history-sync.git
cd codex-history-sync
```

From Gitee:

```bash
git clone git@gitee.com:duke/codex-history-sync.git
cd codex-history-sync
```

## Usage

Check that the CLI starts correctly:

```bash
python3 sync_backend.py --help
```

You can also use the equivalent `help` command:

```bash
python3 sync_backend.py help
```

On Windows PowerShell:

```powershell
py -3 .\sync_backend.py --help
py -3 .\sync_backend.py help
```

Recommended workflow:

### 1. Inspect the Detected Codex Paths

```bash
python3 sync_backend.py status
```

On Windows PowerShell:

```powershell
py -3 .\sync_backend.py status
```

This command does not modify files. It shows:

- The Codex data directory being used
- Whether `config.toml` was found
- Whether `state_5.sqlite` was found
- Whether the `sessions` directory was found
- Whether `session_index.jsonl` was found
- The `model_provider` and `model` that will be synchronized

### 2. Run a Dry Run First

```bash
python3 sync_backend.py sync --dry-run
```

On Windows PowerShell:

```powershell
py -3 .\sync_backend.py sync --dry-run
```

`--dry-run` does not write files. It only reports what would be updated.

### 3. Synchronize History

```bash
python3 sync_backend.py sync
```

On Windows PowerShell:

```powershell
py -3 .\sync_backend.py sync
```

The command creates a backup first, then updates:

- `~/.codex/state_5.sqlite`
- `~/.codex/sessions/**/rollout-*.jsonl`
- `~/.codex/session_index.jsonl`

When synchronizing `session_index.jsonl`, the tool preserves these historical context fields from the database: `cwd`, `git_branch`, `git_sha`, `git_origin_url`, and `rollout_path`. It also fills `git.branch`, `git.commit_hash`, and `git.repository_url`.

Example output:

```text
Synced: provider=openai, model=gpt-5
Database threads: 12/12 updated
Rollout files:    12/12 updated
Index rows:       12/12 updated
Backup: /home/you/.codex/history-sync-backups/codex-history-20260513-230037.tar.gz
```

By default, the tool reads `CODEX_HOME`; if it is not set, it uses `~/.codex`.

## Custom Codex Data Directory

Use `--codex-home` if your Codex data is not in the default location:

```bash
python3 sync_backend.py --codex-home ~/.codex sync
```

On Windows PowerShell:

```powershell
py -3 .\sync_backend.py --codex-home "$env:USERPROFILE\.codex" sync
```

You can also use an environment variable:

```bash
export CODEX_HOME=/path/to/codex-home
python3 sync_backend.py sync
```

On Windows PowerShell:

```powershell
$env:CODEX_HOME = "$env:USERPROFILE\.codex"
py -3 .\sync_backend.py sync
```

## Commands

On Windows PowerShell, replace `python3 sync_backend.py` with `py -3 .\sync_backend.py`.

| Command | Description | Writes Files |
| --- | --- | --- |
| `python3 sync_backend.py help` | Show global help | No |
| `python3 sync_backend.py status` | Show current config and detected paths | No |
| `python3 sync_backend.py sync --dry-run` | Preview what would be synchronized | No |
| `python3 sync_backend.py sync` | Synchronize Codex history | Yes |
| `python3 sync_backend.py backup` | Create a backup manually | Yes, backup only |
| `python3 sync_backend.py list-backups` | List available backups | No |
| `python3 sync_backend.py restore <backup>` | Restore from a backup archive | Yes |

## Backup and Restore

Create a backup:

```bash
python3 sync_backend.py backup
```

List backups:

```bash
python3 sync_backend.py list-backups
```

Restore a backup:

```bash
python3 sync_backend.py restore ~/.codex/history-sync-backups/codex-history-YYYYMMDD-HHMMSS.tar.gz
```

Replace the path with a real backup path from `list-backups`.

## Safety

- `sync --dry-run` only reports changes and does not write files.
- `sync` backs up `config.toml`, `state_5.sqlite`, `session_index.jsonl`, and all `rollout-*.jsonl` files before writing.
- Rebuilding `session_index.jsonl` preserves extra existing fields where possible, including Git branch metadata.
- JSONL files are written atomically to avoid truncation on failure.
- `restore` refuses to extract paths outside `CODEX_HOME`.

## Development and Testing

Run the unit tests:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

On Windows PowerShell:

```powershell
py -3 -m unittest discover -s tests -v
```

Check CLI help:

```bash
python3 sync_backend.py --help
```

## Common Paths

Windows default Codex data directory:

```text
C:\Users\you\.codex
```

Linux/Ubuntu and macOS default Codex data directory:

```text
~/.codex
```

For a custom Codex data directory:

```bash
export CODEX_HOME=/path/to/codex-home
python3 sync_backend.py sync
```

## License

MIT
