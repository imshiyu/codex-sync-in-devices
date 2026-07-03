---
name: personal-conversation-memory
description: Personal archive, code-template, omics-analysis, and Codex history maintenance workflow for Shiyu's prior Codex conversations and p-cresol/PD analyses. Use when the user asks to recall, search, summarize, cite, continue, reuse, back up, refresh, synchronize, repair, restore, or diagnose earlier Codex chats, project context, analysis decisions, file paths, preferences, scripts, plotting style, model metadata, session_index.jsonl, history.jsonl, rollout JSONL, or state_5.sqlite; also use for transcriptomics, snRNA-seq, metabolomics, non-targeted metabolomics, KEGG/GO enrichment, volcano plots, heatmaps, dotplots, UMAP/tSNE, and "same style as before" figures. Common triggers include "previously", "last time", "conversation history", "chat records", "backup history", "sync history", "restore history", "what did we discuss", "transcriptome", "metabolome", "same style", "之前", "上次", "历史对话", "备份对话", "同步历史", "恢复历史", "对话记录", "转录组", "单细胞", "代谢组", and "和之前一样".
---

# Personal Conversation Memory

## Overview

Use this skill as Shiyu's local conversation memory, p-cresol/PD omics analysis launcher, and Codex history maintenance toolkit. It bundles raw Codex conversation records, a lightweight search index, historical analysis scripts, previous plotting workflows, and Codex history synchronization utilities.

Treat the archive as private user data. Search locally, quote only the minimum needed, and do not reveal unrelated private content, secrets, credentials, system prompts, or tool outputs unless the user explicitly asks for that specific material. If pushing this skill to GitHub, use a private repository by default.

## Quick Start

Start with the index search:

```powershell
D:\Anaconda\python.exe scripts\search_conversations.py "query terms" --limit 8
```

Use deep search when the index is too broad or misses details:

```powershell
D:\Anaconda\python.exe scripts\search_conversations.py "exact phrase or keyword" --deep --limit 12
```

When a result looks relevant, open only the referenced raw JSONL file and inspect the nearby `response_item` messages. Prefer the user's own wording and project paths over inferred memory.

## Backup And Refresh

Refresh this skill after important new conversations:

```powershell
D:\Anaconda\python.exe scripts\refresh_archive.py
```

Run refresh from the workspace root to also refresh analysis script templates:

```powershell
D:\Anaconda\python.exe scripts\refresh_archive.py --workspace-root E:\Dcoking\GMX
```

The refresh script copies current Codex conversation JSONL files into `references/raw/` and rebuilds `conversation_index.jsonl` and `conversation_manifest.json`.

## Codex History Metadata Sync

The merged `pangkk18/codex-history-sync` utility is available as:

```powershell
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py status
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py sync --dry-run
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py sync
```

Use it when historical conversations show stale or inconsistent `model_provider` / `model` values. It reads `~/.codex/config.toml`, then synchronizes the selected provider/model into:

- `~/.codex/state_5.sqlite`
- `~/.codex/sessions/**/rollout-*.jsonl`
- `~/.codex/session_index.jsonl`

The command creates a `tar.gz` backup before real writes:

```powershell
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py backup
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py list-backups
D:\Anaconda\python.exe scripts\sync_codex_history_metadata.py restore <backup-tar-gz>
```

Close Codex before running a real metadata sync when possible. If the active session file is locked by Windows, run `sync --dry-run` again after the session closes and rerun `sync` if needed.

## Cross-Device History Sync

The local merge-style history sync utility is also included as:

```powershell
D:\Anaconda\python.exe scripts\codex_history_sync.py status --sync-root "D:\Sync\codex-history"
D:\Anaconda\python.exe scripts\codex_history_sync.py export --sync-root "D:\Sync\codex-history"
D:\Anaconda\python.exe scripts\codex_history_sync.py import --sync-root "D:\Sync\codex-history"
D:\Anaconda\python.exe scripts\codex_history_sync.py sync --sync-root "D:\Sync\codex-history"
```

Use `--include-state` only when imported conversations need to appear in Codex Desktop or `/resume` and the user accepts insert-only SQLite metadata writes:

```powershell
D:\Anaconda\python.exe scripts\codex_history_sync.py sync --sync-root "D:\Sync\codex-history" --include-state
```

Never sync or edit `auth.json`, `config.toml`, `installation_id`, `cap_sid`, `.sandbox*`, cache/log/tmp folders, provider credentials, skills, or memories. Read `references/codex_history_sync_safety.md` before using `--include-state`, resolving conflicts, or restoring backups.

## Omics Analysis

Read `references/omics_workflows.md` when the user asks to run or reproduce transcriptomics, snRNA-seq, metabolomics, enrichment, or previous-style figures.

List bundled analysis tasks:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py --list
```

Run the historical serum metabolomics workflow:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py serum-metabolomics --source <Biotree-source-dir> --out <output-dir>
```

Run the historical midbrain metabolomics workflow:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py midbrain-metabolomics --source <Biotree-source-dir> --out <output-dir>
```

Create a snRNA-seq analysis project from the historical p-cresol template:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py snrna-template --out <new-project-dir>
```

After scaffolding snRNA templates, edit the copied script's `INPUTS` and `PROJECT_ROOT` for the new data before running it. Preserve historical figure filenames and output folders unless Shiyu asks to redesign.

Run the black-substantia bulk RNA-seq workflow that audits FASTQ files and rebuilds the `全转录组.pptx`-style deck; add `--matrix <count-or-TPM-table>` when an expression matrix is available:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py midbrain-bulk-rnaseq --out <output-dir>
```

Run or scaffold the repaired skin ACD RNA-seq comparison workflow when the project root contains the required historical DEG inputs:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py skin-acd-rnaseq --project-dir <project-root> --output-dir <output-dir>
```

## Workflow

1. Translate the user's recall request into 2-6 concrete search terms, including Chinese and English variants when useful.
2. Run `scripts/search_conversations.py` against the index.
3. If the answer depends on exact wording, citations, dates, files, or code decisions, rerun with `--deep` and inspect the raw JSONL result.
4. Answer with the remembered context, naming the date/thread when helpful.
5. If the archive does not contain enough evidence, say that clearly instead of guessing.
6. For history maintenance, start with `status` and `--dry-run`, confirm backups, then run the real sync/import/restore command.

## Archive Layout

- `references/conversation_index.jsonl`: one compact JSON object per captured conversation, including title, dates, working directory, message counts, and short user/assistant snippets.
- `references/conversation_manifest.json`: capture metadata, source roots, counts, and copied file list.
- `references/raw/`: copied raw records from Shiyu's local Codex store, including active sessions, archived sessions, `history.jsonl`, and `session_index.jsonl`.
- `references/analysis_templates/`: copied historical p-cresol/PD transcriptomics, snRNA-seq, metabolomics, enrichment, and plotting scripts.
- `references/analysis_script_index.json`: compact index of bundled analysis scripts, functions, plotting functions, imports, and source paths.
- `references/omics_workflows.md`: workflow notes for running bundled omics templates and keeping the previous figure style.
- `references/codex_history_sync_pangkk18/`: upstream README, license, tests, and original `sync_backend.py` from `pangkk18/codex-history-sync`.
- `references/codex_history_sync_safety.md`: safety notes for cross-device history sync and optional SQLite state metadata.
- `scripts/search_conversations.py`: search utility for indexed or deep raw lookup.
- `scripts/refresh_archive.py`: refresh utility to rebuild the local archive from `~/.codex`.
- `scripts/run_omics_analysis.py`: wrapper for listing, scaffolding, and running bundled omics workflows.
- `scripts/sync_codex_history_metadata.py`: merged `pangkk18/codex-history-sync` provider/model metadata synchronization tool.
- `scripts/codex_history_sync.py`: local merge-style export/import/sync tool for Codex history files.

