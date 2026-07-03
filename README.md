# personal-conversation-memory

Merged private Codex skill for Shiyu's conversation archive, omics workflow reuse, and Codex history maintenance.

This repository contains private conversation history under `references/raw/`. Keep the GitHub repository private unless those files are intentionally removed.

## Main Commands

```powershell
D:\Anaconda\python.exe scripts\search_conversations.py "query terms" --limit 8
D:\Anaconda\python.exe scripts\refresh_archive.py --workspace-root E:\Dcoking\GMX
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py sync --dry-run
D:\Anaconda\python.exe scripts\codex_history_sync.py status --sync-root "D:\Sync\codex-history"
```

See `SKILL.md` for the full operating instructions.

