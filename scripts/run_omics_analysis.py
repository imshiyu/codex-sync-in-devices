#!/usr/bin/env python3
"""Run or scaffold Shiyu's historical omics analysis templates."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
TEMPLATE_DIR = SKILL_DIR / "references" / "analysis_templates"
SCRIPT_INDEX = SKILL_DIR / "references" / "analysis_script_index.json"


TASKS = {
    "serum-metabolomics": {
        "script": "metabolomics/rerun_serum_nontarget_ctrl_pd.py",
        "description": "Run serum non-targeted metabolomics Ctrl vs PD reanalysis.",
        "requires": ["--source", "--out"],
    },
    "midbrain-metabolomics": {
        "script": "metabolomics/rerun_midbrain_nontarget_ctrl_pd.py",
        "description": "Run midbrain non-targeted metabolomics Ctrl vs PD reanalysis.",
        "requires": ["--source", "--out"],
    },
    "filtered-metabolomics-enrichment": {
        "script": "metabolomics/run_filtered_nontarget_enrichment.py",
        "description": "Run the historical post-filtering KEGG/class enrichment workflow.",
        "requires": [],
    },
    "snrna-template": {
        "script": "transcriptomics/p_cresol_snrna_full_analysis.py",
        "description": "Copy the full p-cresol snRNA-seq analysis template for editing input paths.",
        "requires": ["--out"],
        "scaffold_only": True,
    },
    "midbrain-bulk-rnaseq": {
        "script": "transcriptomics/analyze_midbrain_bulk_rnaseq.py",
        "description": "Run or scaffold the black-substantia bulk RNA-seq FASTQ audit, plots, and same-style PPT workflow.",
        "requires": [],
    },
    "skin-acd-rnaseq": {
        "script": "transcriptomics/skin_acd_analysis_repaired.py",
        "description": "Run or scaffold the repaired skin ACD RNA-seq comparison workflow for DMF/OXA/URU Th-polarization analysis.",
        "requires": [],
    },
}


SUPPORT_FILES = {
    "serum-metabolomics": [],
    "midbrain-metabolomics": ["metabolomics/rerun_serum_nontarget_ctrl_pd.py"],
    "filtered-metabolomics-enrichment": [],
    "snrna-template": [
        "transcriptomics/plot_astrocyte_b2m_heatmap_violin.py",
        "transcriptomics/plot_astrocyte_refined_marker_dotplot.py",
        "transcriptomics/plot_astrocyte_refined_top_marker_dotplot.py",
        "transcriptomics/run_astrocyte_deg_table_enrichment.py",
    ],
    "midbrain-bulk-rnaseq": [],
    "skin-acd-rnaseq": [],
}


def rel_to_template(path: str) -> Path:
    return TEMPLATE_DIR / path


def copy_task_files(task_name: str, outdir: Path) -> Path:
    task = TASKS[task_name]
    outdir.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for rel in [task["script"], *SUPPORT_FILES.get(task_name, [])]:
        src = rel_to_template(rel)
        dest = outdir / Path(rel).name
        shutil.copy2(src, dest)
        copied.append(dest.name)
    manifest = {
        "task": task_name,
        "description": task["description"],
        "copied": copied,
        "notes": [
            "These are historical templates copied from Shiyu's prior analyses.",
            "Keep output folder naming and figure filenames unless the user asks for a new style.",
        ],
    }
    (outdir / "template_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return outdir / Path(task["script"]).name


def list_tasks() -> None:
    print("Available tasks:")
    for name, task in TASKS.items():
        print(f"- {name}: {task['description']}")
        print(f"  script: {task['script']}")
        if task["requires"]:
            print(f"  requires: {' '.join(task['requires'])}")
    if SCRIPT_INDEX.exists():
        index = json.loads(SCRIPT_INDEX.read_text(encoding="utf-8"))
        print(f"\nIndexed historical scripts: {len(index.get('scripts', []))}")


def run_python(script: Path, extra_args: list[str], cwd: Path | None = None) -> int:
    cmd = [sys.executable, str(script), *extra_args]
    print("Running:", " ".join(cmd), flush=True)
    return subprocess.call(cmd, cwd=str(cwd) if cwd else None)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", nargs="?", choices=sorted(TASKS), help="Analysis task to run or scaffold.")
    parser.add_argument("--list", action="store_true", help="List available tasks.")
    parser.add_argument("--source", type=Path, help="Input analysis source directory.")
    parser.add_argument("--out", type=Path, help="Output directory.")
    parser.add_argument("--scaffold", action="store_true", help="Copy template scripts to --out without running.")
    parser.add_argument("--report-only", action="store_true", help="Pass --report-only to supported metabolomics tasks.")
    args, unknown = parser.parse_known_args()

    if args.list or not args.task:
        list_tasks()
        return 0

    task = TASKS[args.task]
    if "--source" in task["requires"] and not args.source:
        parser.error(f"{args.task} requires --source")
    if "--out" in task["requires"] and not args.out:
        parser.error(f"{args.task} requires --out")

    scaffold_only = args.scaffold or bool(task.get("scaffold_only"))
    if scaffold_only:
        if not args.out:
            parser.error("--out is required when scaffolding templates")
        script = copy_task_files(args.task, args.out)
        print(f"Scaffolded {args.task} to {args.out}")
        print(f"Main script: {script}")
        if task.get("scaffold_only"):
            print("Edit hard-coded input paths in the copied script, then run it in that project folder.")
        return 0

    script = rel_to_template(task["script"])
    extra_args: list[str] = []
    if args.source:
        extra_args.extend(["--source", str(args.source)])
    if args.out:
        extra_args.extend(["--out", str(args.out)])
    if args.report_only:
        extra_args.append("--report-only")
    extra_args.extend(unknown)
    return run_python(script, extra_args)


if __name__ == "__main__":
    raise SystemExit(main())
