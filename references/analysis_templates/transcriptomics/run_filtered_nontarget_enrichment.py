#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import annotations

import math
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.cm as mcm
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import hypergeom


DATASETS = [
    {
        "name": "serum",
        "label": "血清非靶",
        "source": Path(r"E:\project\血清非靶"),
        "reanalysis": Path(r"E:\project\血清非靶_Ctrl_vs_PD_reanalysis"),
        "kegg_sheet": "Ctrl vs ppm30",
        "removed": ["Ctrl1/SC1", "Ctrl6/SC6", "30ppm4/SP4", "30ppm6/SP6"],
    },
    {
        "name": "midbrain",
        "label": "中脑非靶",
        "source": Path(r"E:\project\中脑非靶"),
        "reanalysis": Path(r"E:\project\中脑非靶_Ctrl_vs_PD_reanalysis_drop_Ctrl2_Ctrl3_Ctrl8_PD3_PD4_PD5"),
        "kegg_sheet": "Ctrl vs PD",
        "removed": ["Ctrl2", "Ctrl3", "Ctrl8", "PD3", "PD4", "PD5"],
    },
]

HEATMAP_COLORS = [
    "#D73027",
    "#F46D43",
    "#FDAE61",
    "#FEE090",
    "#E0F3F8",
    "#ABD9E9",
    "#74ADD1",
    "#4575B4",
]
DOT_CMAP = LinearSegmentedColormap.from_list("filtered_enrichment", HEATMAP_COLORS)

CLASS_LEVELS = ["Super.Class", "Class", "Sub.Class"]

AVAILABLE_FONTS = {font.name for font in fm.fontManager.ttflist}
for font_name in ["Microsoft YaHei", "SimHei", "SimSun", "DejaVu Sans"]:
    if font_name in AVAILABLE_FONTS:
        plt.rcParams["font.family"] = font_name
        break
plt.rcParams["axes.unicode_minus"] = False


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_sheet_name(name: str) -> str:
    return re.sub(r"[\[\]\:\*\?\/\\]", "_", str(name))[:31] or "Sheet"


def safe_filename(value: object) -> str:
    text = str(value)
    text = re.sub(r"[\\/:*?\"<>|]+", "_", text)
    text = re.sub(r"\s+", "_", text).strip("_")
    return text[:120] or "plot"


def wrap_label(value: object, width: int = 48) -> str:
    text = str(value)
    text = re.sub(r"\s+-\s+Mus musculus.*$", "", text)
    text = text.replace("_", " ")
    if text.isupper():
        text = text.title()
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def split_kegg_ids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item for item in re.findall(r"C\d{5}", str(value)) if item]


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    q[valid] = out
    return q


def feature_label(row: pd.Series) -> str:
    for col in ["MS2 name", "MS1 name"]:
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    if pd.notna(row.get("ID")):
        return f"feature_{row['ID']}"
    return "feature"


def prepare_stat_table(stat_file: Path) -> pd.DataFrame:
    df = pd.read_excel(stat_file, sheet_name="Ctrl vs PD")
    df["Diff"] = df["Diff"].astype(str).str.lower().isin(["true", "1", "yes"])
    for col in ["P-Value", "Q-Value", "VIP", "Fold_Change", "Log_Foldchange"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["Direction"] = "All_DEGs"
    df.loc[df["Diff"] & (df["Log_Foldchange"] > 0), "Direction"] = "Up"
    df.loc[df["Diff"] & (df["Log_Foldchange"] < 0), "Direction"] = "Down"
    df["Feature_Label"] = df.apply(feature_label, axis=1)
    return df


def build_kegg_feature_maps(df: pd.DataFrame) -> tuple[set[str], dict[str, list[str]], set[str], set[str], set[str]]:
    background_ids: set[str] = set()
    id_to_features: dict[str, list[str]] = {}
    up_ids: set[str] = set()
    down_ids: set[str] = set()
    all_diff_ids: set[str] = set()
    for _, row in df.iterrows():
        ids = split_kegg_ids(row.get("KEGG COMPOUND ID"))
        label = feature_label(row)
        for kid in ids:
            background_ids.add(kid)
            id_to_features.setdefault(kid, []).append(label)
            if row.get("Diff") is True or bool(row.get("Diff")):
                all_diff_ids.add(kid)
                if row.get("Log_Foldchange", np.nan) > 0:
                    up_ids.add(kid)
                elif row.get("Log_Foldchange", np.nan) < 0:
                    down_ids.add(kid)
    return background_ids, id_to_features, all_diff_ids, up_ids, down_ids


def run_kegg_enrichment(df: pd.DataFrame, source: Path, sheet_name: str) -> pd.DataFrame:
    path_file = source / "09.KEGG Analysis" / "KEGG Pathway.xlsx"
    if not path_file.exists():
        return pd.DataFrame()
    path_df = pd.read_excel(path_file, sheet_name=sheet_name)
    background_ids, id_to_features, all_diff_ids, up_ids, down_ids = build_kegg_feature_maps(df)
    groups = {
        "All_DEGs": all_diff_ids,
        "Up": up_ids,
        "Down": down_ids,
    }
    rows = []
    m = len(background_ids)
    for group, query_ids in groups.items():
        n = len(query_ids)
        if m == 0 or n == 0:
            continue
        for _, row in path_df.iterrows():
            pathway_ids = set(split_kegg_ids(row.get("Compounds (all)", "")))
            if not pathway_ids:
                pathway_ids = set(split_kegg_ids(row.get("Compounds.(all)", "")))
            pathway_bg = pathway_ids & background_ids
            overlap = pathway_bg & query_ids
            if not overlap:
                continue
            x = len(overlap)
            k = len(pathway_bg)
            p_value = float(hypergeom.sf(x - 1, m, k, n)) if k > 0 else np.nan
            overlap_up = overlap & up_ids
            overlap_down = overlap & down_ids
            overlap_features = []
            for kid in sorted(overlap):
                overlap_features.extend(id_to_features.get(kid, []))
            rows.append(
                {
                    "QueryGroup": group,
                    "Pathway": row.get("Pathway"),
                    "Description": str(row.get("Description", "")).replace(" - Mus musculus (house mouse)", ""),
                    "Class": row.get("Class", ""),
                    "Background_compounds": k,
                    "Query_compounds": n,
                    "Overlap_count": x,
                    "Overlap": f"{x}/{k}",
                    "Compounds_all": ";".join(sorted(pathway_bg)),
                    "Compounds_overlap": ";".join(sorted(overlap)),
                    "Overlap_features": "; ".join(sorted(set(overlap_features))[:80]),
                    "Rich_factor": x / k if k else np.nan,
                    "GeneRatio": x / n if n else np.nan,
                    "P-value": p_value,
                    "Adjusted P-value": p_value,
                    "Up_count": len(overlap_up),
                    "Down_count": len(overlap_down),
                    "DA_score": (len(overlap_up) - len(overlap_down)) / x if x else np.nan,
                }
            )
    enrich = pd.DataFrame(rows)
    if enrich.empty:
        return enrich
    frames = []
    for group, group_df in enrich.groupby("QueryGroup", sort=False):
        group_df = group_df.copy()
        group_df["Adjusted P-value"] = bh_fdr(group_df["P-value"].to_numpy())
        frames.append(group_df)
    enrich = pd.concat(frames, ignore_index=True)
    enrich = enrich.sort_values(["QueryGroup", "Adjusted P-value", "P-value", "Overlap_count"], ascending=[True, True, True, False])
    return enrich


def run_class_enrichment(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_masks = {
        "All_DEGs": df["Diff"],
        "Up": df["Diff"] & (df["Log_Foldchange"] > 0),
        "Down": df["Diff"] & (df["Log_Foldchange"] < 0),
    }
    for level in CLASS_LEVELS:
        if level not in df.columns:
            continue
        background = df[df[level].notna() & (df[level].astype(str).str.strip() != "")].copy()
        if background.empty:
            continue
        m = len(background)
        for group, mask in group_masks.items():
            query = df[mask & df[level].notna() & (df[level].astype(str).str.strip() != "")].copy()
            n = len(query)
            if n == 0:
                continue
            for term, term_bg in background.groupby(level, dropna=True):
                term = str(term).strip()
                if not term or term.lower() == "nan":
                    continue
                k = len(term_bg)
                term_query = query[query[level].astype(str) == term]
                x = len(term_query)
                if x == 0:
                    continue
                p_value = float(hypergeom.sf(x - 1, m, k, n))
                labels = term_query["Feature_Label"].dropna().astype(str).drop_duplicates().head(80).tolist()
                rows.append(
                    {
                        "QueryGroup": group,
                        "ClassLevel": level,
                        "Term": term,
                        "Background_features": k,
                        "Query_features": n,
                        "Overlap_count": x,
                        "Overlap": f"{x}/{k}",
                        "Rich_factor": x / k if k else np.nan,
                        "GeneRatio": x / n if n else np.nan,
                        "Fold_enrichment": (x / n) / (k / m) if k and n and m else np.nan,
                        "P-value": p_value,
                        "Adjusted P-value": p_value,
                        "Overlap_features": "; ".join(labels),
                    }
                )
    enrich = pd.DataFrame(rows)
    if enrich.empty:
        return enrich
    frames = []
    for (group, level), group_df in enrich.groupby(["QueryGroup", "ClassLevel"], sort=False):
        group_df = group_df.copy()
        group_df["Adjusted P-value"] = bh_fdr(group_df["P-value"].to_numpy())
        frames.append(group_df)
    enrich = pd.concat(frames, ignore_index=True)
    enrich = enrich.sort_values(["ClassLevel", "QueryGroup", "Adjusted P-value", "P-value", "Overlap_count"], ascending=[True, True, True, True, False])
    return enrich


def plot_enrichment_dotplot(df: pd.DataFrame, out_prefix: Path, title: str, label_col: str = "Description", top_n: int = 20) -> None:
    if df.empty:
        return
    plot_df = df.copy()
    sig = plot_df[plot_df["Adjusted P-value"] <= 0.05].copy()
    plot_df = sig.head(top_n).copy() if len(sig) >= 5 else plot_df.head(top_n).copy()
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values(["Rich_factor", "Overlap_count"], ascending=[True, True]).reset_index(drop=True)
    labels = [wrap_label(value) for value in plot_df[label_col]]
    height = max(5.6, len(plot_df) * 0.52)
    fig, ax = plt.subplots(figsize=(10.8, height))
    fig.subplots_adjust(left=0.44, right=0.86, top=0.91, bottom=0.12)
    p_adjust = plot_df["Adjusted P-value"].astype(float).replace([np.inf, -np.inf], np.nan).fillna(1.0)
    vmin = float(p_adjust.min())
    vmax = float(p_adjust.max())
    if np.isclose(vmin, vmax):
        vmax = vmin + 1e-6
    norm = Normalize(vmin=vmin, vmax=vmax)
    colors = [DOT_CMAP(norm(value)) for value in p_adjust]
    count_min = int(plot_df["Overlap_count"].min())
    count_max = int(plot_df["Overlap_count"].max())
    if count_min == count_max:
        sizes = np.array([220] * len(plot_df))
    else:
        sizes = 65 + ((plot_df["Overlap_count"] - count_min) / float(count_max - count_min)) * 380
    ax.scatter(plot_df["Rich_factor"], range(len(plot_df)), s=sizes, c=colors, edgecolors="#666666", linewidth=0.5, zorder=4)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Rich factor")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.grid(True, color="#CCCCCC", linewidth=0.5, zorder=0)
    for spine in ax.spines.values():
        spine.set_color("#666666")
        spine.set_linewidth(0.8)
    scalar = mcm.ScalarMappable(cmap=DOT_CMAP, norm=norm)
    scalar.set_array([])
    cbar = plt.colorbar(scalar, ax=ax, shrink=0.35, pad=0.03, aspect=12)
    cbar.set_label("p.adjust", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    for suffix in [".png", ".pdf"]:
        fig.savefig(str(out_prefix) + suffix, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_da_score(kegg_df: pd.DataFrame, out_prefix: Path, title: str, top_n: int = 20) -> None:
    if kegg_df.empty:
        return
    plot_df = kegg_df[kegg_df["QueryGroup"] == "All_DEGs"].copy()
    if plot_df.empty:
        return
    sig = plot_df[plot_df["Adjusted P-value"] <= 0.05].copy()
    plot_df = sig.head(top_n).copy() if len(sig) >= 5 else plot_df.head(top_n).copy()
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values("DA_score", ascending=True)
    labels = [wrap_label(value, width=50) for value in plot_df["Description"]]
    colors = np.where(plot_df["DA_score"] >= 0, "#D73027", "#4575B4")
    height = max(5.6, len(plot_df) * 0.5)
    fig, ax = plt.subplots(figsize=(10.5, height))
    fig.subplots_adjust(left=0.45, right=0.96, top=0.91, bottom=0.12)
    ax.barh(range(len(plot_df)), plot_df["DA_score"], color=colors, edgecolor="#666666", linewidth=0.4)
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("DA score = (up - down) / overlap")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.grid(axis="x", color="#CCCCCC", linewidth=0.5)
    for suffix in [".png", ".pdf"]:
        fig.savefig(str(out_prefix) + suffix, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_class_barplot(class_df: pd.DataFrame, out_prefix: Path, title: str, group: str, level: str, top_n: int = 18) -> None:
    sub = class_df[(class_df["QueryGroup"] == group) & (class_df["ClassLevel"] == level)].copy()
    if sub.empty:
        return
    sig = sub[sub["Adjusted P-value"] <= 0.05].copy()
    sub = sig.head(top_n).copy() if len(sig) >= 5 else sub.head(top_n).copy()
    if sub.empty:
        return
    sub["minus_log10_padj"] = -np.log10(sub["Adjusted P-value"].clip(lower=1e-300))
    sub = sub.sort_values("minus_log10_padj", ascending=True)
    labels = [wrap_label(value, width=44) for value in sub["Term"]]
    height = max(5.4, len(sub) * 0.48)
    fig, ax = plt.subplots(figsize=(9.5, height))
    fig.subplots_adjust(left=0.43, right=0.96, top=0.91, bottom=0.12)
    ax.barh(range(len(sub)), sub["minus_log10_padj"], color="#74ADD1", edgecolor="#666666", linewidth=0.45)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("-log10(adjusted P-value)")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.grid(axis="x", color="#CCCCCC", linewidth=0.5)
    for suffix in [".png", ".pdf"]:
        fig.savefig(str(out_prefix) + suffix, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_excel(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=safe_sheet_name(name), index=False)


def write_summary(dataset: dict, outdir: Path, stat_df: pd.DataFrame, kegg: pd.DataFrame, class_enrich: pd.DataFrame) -> None:
    diff = stat_df[stat_df["Diff"]]
    rows = [
        ["Total features", len(stat_df)],
        ["Differential features", len(diff)],
        ["PD up features", int((diff["Log_Foldchange"] > 0).sum())],
        ["PD down features", int((diff["Log_Foldchange"] < 0).sum())],
        ["KEGG enrichment rows", len(kegg)],
        ["KEGG FDR<=0.05 rows", int((kegg["Adjusted P-value"] <= 0.05).sum()) if not kegg.empty else 0],
        ["Class enrichment rows", len(class_enrich)],
        ["Class FDR<=0.05 rows", int((class_enrich["Adjusted P-value"] <= 0.05).sum()) if not class_enrich.empty else 0],
    ]
    lines = [f"# {dataset['label']} 删样后富集分析", ""]
    lines.append("- 重分析目录：`{0}`".format(dataset["reanalysis"]))
    lines.append("- 删除样本：{0}".format("、".join(dataset["removed"])))
    lines.append("- 差异筛选口径：使用重分析统计表中的 `Diff=True`，即 VIP>1 且 P-Value<0.05。")
    lines.append("")
    lines.append("| 指标 | 数量 |")
    lines.append("|---|---:|")
    for key, value in rows:
        lines.append(f"| {key} | {value} |")
    lines.append("")
    lines.append("## 输出文件")
    lines.append("")
    lines.append("- `KEGG_pathway_enrichment_extended.xlsx`：All_DEGs、Up、Down 的 KEGG ORA 富集。")
    lines.append("- `Chemical_class_enrichment.xlsx`：Super.Class、Class、Sub.Class 代谢物类别富集。")
    lines.append("- `*_bubble.png/pdf`：KEGG 富集气泡图。")
    lines.append("- `KEGG_DA_score.png/pdf`：通路上调/下调方向得分图。")
    lines.append("- `Class_*_barplot.png/pdf`：代谢物类别富集柱状图。")
    if not kegg.empty:
        lines.append("")
        lines.append("## KEGG Top 10")
        top = kegg.head(10)[["QueryGroup", "Pathway", "Description", "Overlap", "Rich_factor", "P-value", "Adjusted P-value"]].copy()
        lines.extend(markdown_table(top))
    if not class_enrich.empty:
        lines.append("")
        lines.append("## Chemical Class Top 10")
        top = class_enrich.head(10)[["QueryGroup", "ClassLevel", "Term", "Overlap", "Fold_enrichment", "P-value", "Adjusted P-value"]].copy()
        lines.extend(markdown_table(top))
    (outdir / "README_enrichment_after_sample_filtering.md").write_text("\n".join(lines), encoding="utf-8")


def markdown_table(df: pd.DataFrame) -> list[str]:
    if df.empty:
        return ["无。"]
    out = df.copy()
    for col in out.columns:
        if np.issubdtype(out[col].dtype, np.number):
            out[col] = out[col].map(lambda x: f"{x:.4g}" if pd.notna(x) else "")
    out = out.fillna("").astype(str)
    lines = [
        "| " + " | ".join(out.columns) + " |",
        "| " + " | ".join(["---"] * len(out.columns)) + " |",
    ]
    for row in out.values.tolist():
        lines.append("| " + " | ".join(str(cell).replace("|", "/") for cell in row) + " |")
    return lines


def analyze_dataset(dataset: dict) -> dict:
    reanalysis = dataset["reanalysis"]
    stat_file = reanalysis / "02.Statistical Analysis" / "Statistical Analysis Results.xlsx"
    outdir = reanalysis / "10.Enrichment Analysis" / "Ctrl vs PD" / "extended_enrichment"
    ensure_dir(outdir)
    print(f"Running enrichment: {dataset['label']}", flush=True)
    stat_df = prepare_stat_table(stat_file)
    kegg = run_kegg_enrichment(stat_df, dataset["source"], dataset["kegg_sheet"])
    class_enrich = run_class_enrichment(stat_df)

    if not kegg.empty:
        write_excel(
            outdir / "KEGG_pathway_enrichment_extended.xlsx",
            {
                "All_results": kegg,
                "FDR_0.05": kegg[kegg["Adjusted P-value"] <= 0.05],
                "All_DEGs": kegg[kegg["QueryGroup"] == "All_DEGs"],
                "Up": kegg[kegg["QueryGroup"] == "Up"],
                "Down": kegg[kegg["QueryGroup"] == "Down"],
            },
        )
        kegg.to_csv(outdir / "KEGG_pathway_enrichment_extended.csv", index=False, encoding="utf-8-sig")
        for group in ["All_DEGs", "Up", "Down"]:
            sub = kegg[kegg["QueryGroup"] == group].copy()
            plot_enrichment_dotplot(sub, outdir / f"KEGG_{group}_bubble", f"{dataset['label']} {group}: KEGG pathway enrichment")
        plot_da_score(kegg, outdir / "KEGG_DA_score", f"{dataset['label']} KEGG DA score after sample filtering")

    if not class_enrich.empty:
        sheets = {"All_results": class_enrich, "FDR_0.05": class_enrich[class_enrich["Adjusted P-value"] <= 0.05]}
        for level in CLASS_LEVELS:
            sheets[level] = class_enrich[class_enrich["ClassLevel"] == level]
        write_excel(outdir / "Chemical_class_enrichment.xlsx", sheets)
        class_enrich.to_csv(outdir / "Chemical_class_enrichment.csv", index=False, encoding="utf-8-sig")
        for group in ["All_DEGs", "Up", "Down"]:
            for level in CLASS_LEVELS:
                if class_enrich[(class_enrich["QueryGroup"] == group) & (class_enrich["ClassLevel"] == level)].empty:
                    continue
                plot_class_barplot(
                    class_enrich,
                    outdir / f"Class_{group}_{safe_filename(level)}_barplot",
                    f"{dataset['label']} {group}: {level} enrichment",
                    group,
                    level,
                )

    write_summary(dataset, outdir, stat_df, kegg, class_enrich)
    return {
        "Dataset": dataset["label"],
        "Output": str(outdir),
        "Total_features": len(stat_df),
        "Diff_features": int(stat_df["Diff"].sum()),
        "Up_features": int((stat_df["Diff"] & (stat_df["Log_Foldchange"] > 0)).sum()),
        "Down_features": int((stat_df["Diff"] & (stat_df["Log_Foldchange"] < 0)).sum()),
        "KEGG_rows": len(kegg),
        "KEGG_FDR_0.05": int((kegg["Adjusted P-value"] <= 0.05).sum()) if not kegg.empty else 0,
        "Class_rows": len(class_enrich),
        "Class_FDR_0.05": int((class_enrich["Adjusted P-value"] <= 0.05).sum()) if not class_enrich.empty else 0,
    }


def main() -> None:
    summaries = [analyze_dataset(dataset) for dataset in DATASETS]
    summary_df = pd.DataFrame(summaries)
    summary_path = Path(r"E:\project\nontarget_enrichment_after_sample_filtering_summary.xlsx")
    with pd.ExcelWriter(summary_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary", index=False)
    summary_df.to_csv(summary_path.with_suffix(".csv"), index=False, encoding="utf-8-sig")
    print(summary_df.to_string(index=False))
    print("Summary:", summary_path)


if __name__ == "__main__":
    main()
