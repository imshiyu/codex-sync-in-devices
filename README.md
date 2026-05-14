# Codex History Sync

[中文](#中文) | [English](#english)

## 中文

用于在多台设备之间同步 Codex 对话历史的 Codex skill。它通过一个私有共享文件夹中转历史数据，例如 OneDrive、Dropbox、Syncthing、iCloud Drive、私有 Git 仓库或移动硬盘。

这个 skill 会同步 Codex 对话历史和可选的线程可见性 metadata，但会刻意排除登录、账号、设备、provider 和配置状态。

### 安装

把整个 `codex-history-sync` 文件夹复制到 Codex skills 目录。

Windows:

```powershell
Copy-Item -Recurse -Force . "$env:USERPROFILE\.codex\skills\codex-history-sync"
```

macOS/Linux:

```bash
cp -R . ~/.codex/skills/codex-history-sync
```

然后重启 Codex，或打开一个新的 Codex 会话，让 skill 列表重新加载。

### 快速使用

先选择一个私有共享文件夹，例如 `D:\Sync\codex-history`。

查看状态：

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" status --sync-root "D:\Sync\codex-history"
```

在第一台设备导出历史：

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" export --sync-root "D:\Sync\codex-history"
```

在另一台设备导入历史：

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" import --sync-root "D:\Sync\codex-history"
```

日常双向同步：

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history"
```

如果导入后 Codex 历史列表或 `/resume` 里看不到会话，可以同步线程 metadata：

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history" --include-state
```

第一次导入前建议先 dry run：

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history" --dry-run
```

### 同步边界

默认同步：

- `sessions/**`
- `archived_sessions/**`
- `session_index.jsonl`
- `history.jsonl`

使用 `--include-state` 时可选同步：

- 从 `state_5.sqlite` 中导出的、插入式的线程 metadata

永不同步：

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

导入时会自动备份到：

```text
~/.codex/backups_state/history-sync/
```

### 工作方式

共享文件夹只是中转站。Codex 不会直接从共享文件夹读取历史；脚本会把一台设备的 `~/.codex` 历史导出到共享文件夹，再从共享文件夹导入到另一台设备自己的 `~/.codex`。

```text
设备 A 的 ~/.codex  ->  共享文件夹  ->  设备 B 的 ~/.codex
```

### 致谢

受 [Dailin521/codex-provider-sync](https://github.com/Dailin521/codex-provider-sync) 的同步目标和安全边界启发。

## English

A Codex skill for synchronizing Codex conversation history across devices through a private shared folder, such as OneDrive, Dropbox, Syncthing, iCloud Drive, a private Git repository, or an external drive.

The skill syncs Codex conversation history and optional thread visibility metadata. It intentionally excludes login, account, device, provider, and configuration state.

### Install

Copy the entire `codex-history-sync` folder to your Codex skills directory.

Windows:

```powershell
Copy-Item -Recurse -Force . "$env:USERPROFILE\.codex\skills\codex-history-sync"
```

macOS/Linux:

```bash
cp -R . ~/.codex/skills/codex-history-sync
```

Then restart Codex or open a new Codex session so the skill list is refreshed.

### Quick Usage

Pick a private shared folder, for example `D:\Sync\codex-history`.

Check status:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" status --sync-root "D:\Sync\codex-history"
```

Export history from the first device:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" export --sync-root "D:\Sync\codex-history"
```

Import history on another device:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" import --sync-root "D:\Sync\codex-history"
```

Run bidirectional sync for regular use:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history"
```

If imported sessions do not appear in Codex history or `/resume`, include thread metadata:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history" --include-state
```

Before the first import, use dry run:

```powershell
python "$env:USERPROFILE\.codex\skills\codex-history-sync\scripts\codex_history_sync.py" sync --sync-root "D:\Sync\codex-history" --dry-run
```

### Safety Boundary

Synced by default:

- `sessions/**`
- `archived_sessions/**`
- `session_index.jsonl`
- `history.jsonl`

Optional with `--include-state`:

- insert-only selected thread metadata exported from `state_5.sqlite`

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

### How It Works

The shared folder is only a transfer store. Codex does not read history directly from the shared folder; the script exports history from one device's `~/.codex` into the shared folder, then imports it into another device's own `~/.codex`.

```text
Device A ~/.codex  ->  shared folder  ->  Device B ~/.codex
```

### Credits

Inspired by the safety boundary and sync goals of [Dailin521/codex-provider-sync](https://github.com/Dailin521/codex-provider-sync).
