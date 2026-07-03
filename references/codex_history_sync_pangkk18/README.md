# Codex History Sync

[English](README_EN.md) | 简体中文

一个面向 Windows、Linux、Ubuntu 和 macOS 的 Codex 历史会话同步工具。它没有 UI，只提供命令行，适合在服务器、开发机或脚本环境中直接运行。

这个工具用于把当前 Codex 配置中的模型信息同步到历史会话元数据中，解决历史会话列表、索引或会话文件里显示的 `model_provider` / `model` 与当前配置不一致的问题。

Codex 的历史元数据通常散落在这些位置：

- `~/.codex/config.toml`
- `~/.codex/state_5.sqlite`
- `~/.codex/sessions/**/rollout-*.jsonl`
- `~/.codex/session_index.jsonl`

本工具会读取 `config.toml` 中当前的 `model_provider` / `model`，并同步到历史数据库、历史会话 JSONL 和会话索引。每次实际写入前都会自动创建 `tar.gz` 备份。

## 特性

- 支持 Windows、Linux、Ubuntu 和 macOS。
- 无第三方依赖，只需要 Python。
- 自动读取 `CODEX_HOME`，未设置时使用 `~/.codex`。
- 支持同步 SQLite 历史数据库、会话 JSONL 和全局会话索引。
- 重建会话索引时会保留工作目录、Git 分支、提交哈希和远端地址等历史上下文。
- 支持 `--dry-run` 试运行。
- 正式写入前自动创建 `tar.gz` 备份。
- 支持列出备份和从备份恢复。

## 环境要求

- Python 3.9+
- 无第三方依赖

## 获取代码

```bash
git clone git@gitee.com:duke/codex-history-sync.git
cd codex-history-sync
```

## 如何使用

先进入本项目目录：

```bash
cd codex-history-sync
```

查看工具是否可以正常启动：

```bash
python3 sync_backend.py --help
```

也可以使用等价的 `help` 命令：

```bash
python3 sync_backend.py help
```

Windows PowerShell 中也可以使用：

```powershell
py -3 .\sync_backend.py --help
py -3 .\sync_backend.py help
```

推荐按下面顺序使用。

### 1. 查看当前 Codex 路径和模型配置

```bash
python3 sync_backend.py status
```

Windows PowerShell：

```powershell
py -3 .\sync_backend.py status
```

这个命令不会修改任何文件，只会显示：

- 当前使用的 Codex 数据目录
- 是否找到了 `config.toml`
- 是否找到了 `state_5.sqlite`
- 是否找到了 `sessions` 目录
- 是否找到了 `session_index.jsonl`
- 将要同步的 `model_provider` 和 `model`

### 2. 先试运行

```bash
python3 sync_backend.py sync --dry-run
```

Windows PowerShell：

```powershell
py -3 .\sync_backend.py sync --dry-run
```

`--dry-run` 不会写入任何文件，只会统计哪些内容需要更新。建议第一次使用时先运行这个命令。

### 3. 正式同步历史会话

```bash
python3 sync_backend.py sync
```

Windows PowerShell：

```powershell
py -3 .\sync_backend.py sync
```

正式同步时，工具会先自动创建备份，然后更新：

- `~/.codex/state_5.sqlite`
- `~/.codex/sessions/**/rollout-*.jsonl`
- `~/.codex/session_index.jsonl`

同步 `session_index.jsonl` 时会从数据库保留这些历史上下文字段：`cwd`、`git_branch`、`git_sha`、`git_origin_url`、`rollout_path`，并补齐 `git.branch`、`git.commit_hash`、`git.repository_url`。

命令输出示例：

```text
Synced: provider=openai, model=gpt-5
Database threads: 12/12 updated
Rollout files:    12/12 updated
Index rows:       12/12 updated
Backup: /home/you/.codex/history-sync-backups/codex-history-20260513-230037.tar.gz
```

默认读取 `CODEX_HOME`，如果没有设置则使用 `~/.codex`。

## 指定 Codex 数据目录

如果你的 Codex 数据不在默认的 `~/.codex`，可以通过 `--codex-home` 指定：

```bash
python3 sync_backend.py --codex-home ~/.codex sync
```

Windows PowerShell：

```powershell
py -3 .\sync_backend.py --codex-home "$env:USERPROFILE\.codex" sync
```

也可以使用环境变量：

```bash
export CODEX_HOME=/path/to/codex-home
python3 sync_backend.py sync
```

Windows PowerShell：

```powershell
$env:CODEX_HOME = "$env:USERPROFILE\.codex"
py -3 .\sync_backend.py sync
```

## 命令说明

Windows PowerShell 中可将下面命令里的 `python3 sync_backend.py` 替换为 `py -3 .\sync_backend.py`。

| 命令 | 作用 | 是否写入文件 |
| --- | --- | --- |
| `python3 sync_backend.py help` | 查看全局帮助 | 否 |
| `python3 sync_backend.py status` | 查看当前配置和文件路径 | 否 |
| `python3 sync_backend.py sync --dry-run` | 试运行同步，统计会修改的内容 | 否 |
| `python3 sync_backend.py sync` | 正式同步历史会话 | 是 |
| `python3 sync_backend.py backup` | 手动创建备份 | 是，只写入备份文件 |
| `python3 sync_backend.py list-backups` | 列出已有备份 | 否 |
| `python3 sync_backend.py restore <backup>` | 从备份恢复 | 是 |

## 备份与恢复

创建备份：

```bash
python3 sync_backend.py backup
```

列出备份：

```bash
python3 sync_backend.py list-backups
```

恢复备份：

```bash
python3 sync_backend.py restore ~/.codex/history-sync-backups/codex-history-YYYYMMDD-HHMMSS.tar.gz
```

恢复时请把路径替换成 `list-backups` 输出的真实备份路径。

## 推荐使用流程

```bash
cd codex-history-sync
python3 sync_backend.py status
python3 sync_backend.py sync --dry-run
python3 sync_backend.py sync
```

如果同步后发现历史会话显示异常，可以用下面流程恢复：

```bash
python3 sync_backend.py list-backups
python3 sync_backend.py restore ~/.codex/history-sync-backups/codex-history-YYYYMMDD-HHMMSS.tar.gz
```

## 安全性

- `sync --dry-run` 只统计会修改的内容，不写入文件。
- `sync` 写入前会备份 `config.toml`、`state_5.sqlite`、`session_index.jsonl` 和所有 `rollout-*.jsonl`。
- `session_index.jsonl` 重建时会尽量保留已有额外字段，包括 Git 分支相关信息。
- JSONL 文件使用原子写入，避免中途失败导致文件被截断。
- `restore` 会拒绝释放备份中指向 `CODEX_HOME` 之外的路径。

## 开发与测试

运行单元测试：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

Windows PowerShell：

```powershell
py -3 -m unittest discover -s tests -v
```

查看命令行帮助：

```bash
python3 sync_backend.py --help
```

## 常见路径

Windows 的默认 Codex 数据目录通常是：

```text
C:\Users\you\.codex
```

Linux/Ubuntu 和 macOS 的默认 Codex 数据目录通常都是：

```text
~/.codex
```

如果你的 Codex 使用了自定义目录，设置环境变量即可：

```bash
export CODEX_HOME=/path/to/codex-home
python3 sync_backend.py sync
```

## License

MIT
