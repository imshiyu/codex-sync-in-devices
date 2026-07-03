# Safety Notes

## What is synced

- `sessions/**`
- `archived_sessions/**`
- `session_index.jsonl`
- `history.jsonl`
- Optional insert-only rows from `state_5.sqlite` tables related to thread visibility.

## What is intentionally excluded

- Login and account files: `auth.json`, `cap_sid`
- Device identity and config: `installation_id`, `config.toml`
- Runtime or cache data: `.sandbox*`, `cache`, `log`, `logs_*.sqlite`, `tmp`
- Skills, rules, memories, and user customizations

## State metadata

`state_5.sqlite` can affect whether imported conversations appear in Codex Desktop or `/resume`. The bundled script exports selected thread metadata and imports it insert-only by default. This avoids replacing a device's local state database while still making missing history rows visible.

Use `--state-conflict use-newer` only when the same thread exists on two devices and the user wants the newer metadata row to win.

## Provider/account caveat

The reference project `codex-provider-sync` documents an important boundary: provider metadata can be synchronized for visibility, but encrypted content may remain bound to the original provider/account. Do not promise that cross-account or cross-provider imported sessions can always be resumed.

## Restore behavior

The script prints a backup directory for import/state operations. Restore from that directory if the user wants to undo an import. Close Codex before restoring state metadata.
