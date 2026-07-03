#!/usr/bin/env python3
"""Build the personal conversation archive bundled with this skill."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
from collections import Counter
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
REFERENCES_DIR = SKILL_DIR / "references"
RAW_DIR = REFERENCES_DIR / "raw"
ANALYSIS_TEMPLATES_DIR = REFERENCES_DIR / "analysis_templates"
UUID_RE = re.compile(
    r"([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})",
    re.IGNORECASE,
)

ANALYSIS_TEMPLATE_FILES = {
    "metabolomics": [
        "rerun_serum_nontarget_ctrl_pd.py",
        "rerun_midbrain_nontarget_ctrl_pd.py",
        "run_filtered_nontarget_enrichment.py",
        "redraw_serum_kegg_enrichment_count_legend.py",
    ],
    "transcriptomics": [
        "analyze_midbrain_bulk_rnaseq.py",
        "p_cresol_snrna_full_analysis.py",
        "p_cresol_snrna_replot_broad_neuron.py",
        "plot_astrocyte_b2m_heatmap_violin.py",
        "plot_astrocyte_refined_marker_dotplot.py",
        "plot_astrocyte_refined_top_marker_dotplot.py",
        "plot_astrocyte_refined_c6_c12_merged_embeddings.py",
        "plot_astrocyte_refined_c6_c12_merged_top_marker_b2m.py",
        "run_astrocyte_deg_table_enrichment.py",
        "run_astrocyte_nominal_p_enrichment.py",
        "run_filtered_nontarget_enrichment.py",
        "skin_acd_analysis_repaired.py",
    ],
}


def default_codex_home() -> Path:
    return Path(os.environ.get("CODEX_HOME") or Path.home() / ".codex")


def safe_read_jsonl(path: Path):
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                yield line_number, json.loads(line)
            except json.JSONDecodeError:
                continue


def extract_text(content: Any) -> str:
    parts: list[str] = []

    def visit(value: Any) -> None:
        if isinstance(value, str):
            parts.append(value)
        elif isinstance(value, dict):
            for key in ("text", "input_text", "output_text"):
                item = value.get(key)
                if isinstance(item, str):
                    parts.append(item)
            nested = value.get("content")
            if nested is not None:
                visit(nested)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(content)
    return normalize_text(" ".join(parts))


def normalize_text(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def snippet(text: str, limit: int = 360) -> str:
    text = normalize_text(text)
    if len(text) <= limit:
        return text
    return text[: limit - 3].rstrip() + "..."


def load_thread_index(codex_home: Path) -> dict[str, dict[str, str]]:
    index_path = codex_home / "session_index.jsonl"
    thread_index: dict[str, dict[str, str]] = {}
    if not index_path.exists():
        return thread_index
    for _, obj in safe_read_jsonl(index_path):
        session_id = obj.get("id")
        if isinstance(session_id, str):
            thread_index[session_id] = {
                "thread_name": str(obj.get("thread_name") or ""),
                "updated_at": str(obj.get("updated_at") or ""),
            }
    return thread_index


def iter_conversation_files(codex_home: Path) -> list[Path]:
    files: list[Path] = []
    for rel_root in ("sessions", "archived_sessions"):
        root = codex_home / rel_root
        if root.exists():
            files.extend(path for path in root.rglob("*.jsonl") if path.is_file())
    for name in ("history.jsonl", "session_index.jsonl"):
        path = codex_home / name
        if path.exists():
            files.append(path)
    return sorted(set(files), key=lambda path: str(path).lower())


def relative_to_codex_home(path: Path, codex_home: Path) -> Path:
    try:
        return path.resolve().relative_to(codex_home.resolve())
    except ValueError:
        return Path("external") / path.name


def parse_rollout(
    source_path: Path,
    raw_rel_path: Path,
    thread_index: dict[str, dict[str, str]],
) -> dict[str, Any]:
    session_id = ""
    created_at = ""
    cwd = ""
    source = ""
    originator = ""
    model_provider = ""
    role_counts: Counter[str] = Counter()
    user_prompts: list[str] = []
    assistant_replies: list[str] = []
    preview_parts: list[str] = []

    match = UUID_RE.search(source_path.name)
    if match:
        session_id = match.group(1)

    for _, obj in safe_read_jsonl(source_path):
        obj_type = obj.get("type")
        payload = obj.get("payload") or {}

        if obj_type == "session_meta":
            session_id = str(payload.get("id") or session_id)
            created_at = str(payload.get("timestamp") or created_at)
            cwd = str(payload.get("cwd") or cwd)
            source = str(payload.get("source") or source)
            originator = str(payload.get("originator") or originator)
            model_provider = str(payload.get("model_provider") or model_provider)
            continue

        if obj_type != "response_item" or payload.get("type") != "message":
            continue

        role = str(payload.get("role") or "unknown")
        text = extract_text(payload.get("content"))
        if not text:
            continue

        role_counts[role] += 1
        if role == "user" and len(user_prompts) < 12:
            user_prompts.append(snippet(text))
        elif role == "assistant" and len(assistant_replies) < 8:
            assistant_replies.append(snippet(text))

        if role in {"user", "assistant"} and len(preview_parts) < 36:
            preview_parts.append(f"{role}: {snippet(text, 500)}")

    thread = thread_index.get(session_id, {})
    stat = source_path.stat()
    updated_at = thread.get("updated_at") or dt.datetime.fromtimestamp(
        stat.st_mtime, tz=dt.timezone.utc
    ).isoformat()

    return {
        "session_id": session_id,
        "thread_name": thread.get("thread_name") or "",
        "created_at": created_at,
        "updated_at": updated_at,
        "cwd": cwd,
        "source": source,
        "originator": originator,
        "model_provider": model_provider,
        "message_counts": dict(role_counts),
        "user_prompts": user_prompts,
        "assistant_replies": assistant_replies,
        "source_file": str(source_path),
        "raw_file": str(raw_rel_path).replace("\\", "/"),
        "bytes": stat.st_size,
        "search_text": snippet(
            " ".join(
                [
                    thread.get("thread_name") or "",
                    cwd,
                    source,
                    originator,
                    *user_prompts,
                    *assistant_replies,
                    *preview_parts,
                ]
            ),
            limit=24000,
        ),
    }


def build_archive(codex_home: Path, skill_dir: Path) -> dict[str, Any]:
    references_dir = skill_dir / "references"
    raw_dir = references_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    thread_index = load_thread_index(codex_home)
    files = iter_conversation_files(codex_home)
    copied_files: list[dict[str, Any]] = []
    index_entries: list[dict[str, Any]] = []

    for source_path in files:
        rel_path = relative_to_codex_home(source_path, codex_home)
        raw_path = raw_dir / rel_path
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, raw_path)

        copied_files.append(
            {
                "source": str(source_path),
                "raw_file": str(raw_path.relative_to(skill_dir)).replace("\\", "/"),
                "bytes": source_path.stat().st_size,
                "modified_at": dt.datetime.fromtimestamp(
                    source_path.stat().st_mtime, tz=dt.timezone.utc
                ).isoformat(),
            }
        )

        if source_path.name.startswith("rollout-") and source_path.suffix == ".jsonl":
            index_entries.append(
                parse_rollout(
                    source_path,
                    raw_path.relative_to(skill_dir),
                    thread_index,
                )
            )

    index_path = references_dir / "conversation_index.jsonl"
    with index_path.open("w", encoding="utf-8", newline="\n") as handle:
        for entry in sorted(index_entries, key=lambda item: item.get("updated_at") or ""):
            handle.write(json.dumps(entry, ensure_ascii=False, sort_keys=True) + "\n")

    manifest = {
        "captured_at": dt.datetime.now(tz=dt.timezone.utc).isoformat(),
        "source_codex_home": str(codex_home),
        "skill_dir": str(skill_dir),
        "archive_root": "references/raw",
        "index_file": "references/conversation_index.jsonl",
        "file_count": len(copied_files),
        "conversation_count": len(index_entries),
        "total_bytes": sum(item["bytes"] for item in copied_files),
        "files": copied_files,
    }

    manifest_path = references_dir / "conversation_manifest.json"
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    return manifest


def script_summary(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    functions = re.findall(r"^def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(", text, flags=re.MULTILINE)
    imports = re.findall(r"^(?:import|from)\s+([A-Za-z0-9_\.]+)", text, flags=re.MULTILINE)
    plot_functions = [name for name in functions if name.startswith("plot_")]
    run_functions = [name for name in functions if name.startswith("run_")]
    has_argparse = "argparse" in text
    has_main = "def main(" in text
    return {
        "name": path.name,
        "bytes": path.stat().st_size,
        "modified_at": dt.datetime.fromtimestamp(path.stat().st_mtime, tz=dt.timezone.utc).isoformat(),
        "functions": functions[:80],
        "plot_functions": plot_functions,
        "run_functions": run_functions,
        "imports": sorted(set(imports))[:80],
        "has_argparse": has_argparse,
        "has_main": has_main,
    }


def build_analysis_templates(workspace_root: Path, skill_dir: Path) -> dict[str, Any]:
    references_dir = skill_dir / "references"
    template_dir = references_dir / "analysis_templates"
    template_dir.mkdir(parents=True, exist_ok=True)

    scripts: list[dict[str, Any]] = []
    copied = 0
    missing: list[str] = []

    for domain, names in ANALYSIS_TEMPLATE_FILES.items():
        domain_dir = template_dir / domain
        domain_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            source = workspace_root / name
            dest = domain_dir / name
            if source.exists():
                shutil.copy2(source, dest)
                summary_source = source
                source_label = str(source)
            elif dest.exists():
                summary_source = dest
                source_label = str(dest)
            else:
                missing.append(str(source))
                continue
            summary = script_summary(summary_source)
            summary.update(
                {
                    "domain": domain,
                    "source": source_label,
                    "template_file": str(dest.relative_to(skill_dir)).replace("\\", "/"),
                }
            )
            scripts.append(summary)
            copied += 1

    index = {
        "captured_at": dt.datetime.now(tz=dt.timezone.utc).isoformat(),
        "workspace_root": str(workspace_root),
        "template_root": "references/analysis_templates",
        "script_count": copied,
        "missing": missing,
        "scripts": sorted(scripts, key=lambda item: (item["domain"], item["name"])),
    }

    index_path = references_dir / "analysis_script_index.json"
    with index_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(index, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return index


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--codex-home",
        type=Path,
        default=default_codex_home(),
        help="Codex home directory containing sessions and archived_sessions.",
    )
    parser.add_argument(
        "--skill-dir",
        type=Path,
        default=SKILL_DIR,
        help="Skill directory to populate.",
    )
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path.cwd(),
        help="Workspace root containing historical analysis scripts.",
    )
    args = parser.parse_args()

    manifest = build_archive(args.codex_home.expanduser(), args.skill_dir.resolve())
    script_index = build_analysis_templates(args.workspace_root.resolve(), args.skill_dir.resolve())
    print(
        "Captured {conversation_count} conversations from {file_count} files "
        "({total_bytes} bytes).".format(**manifest)
    )
    print(f"Index: {manifest['index_file']}")
    print(f"Analysis templates: {script_index['script_count']} scripts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
