from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom
from sklearn.decomposition import PCA

import rerun_serum_nontarget_ctrl_pd as base


RAW_CTRL = [f"Ctrl{i}" for i in range(1, 9)]
RAW_PD = [f"PD{i}" for i in range(1, 9)]
DROP_RAW = ["Ctrl2", "Ctrl3", "Ctrl8", "PD3", "PD4", "PD5"]
KEEP_CTRL = ["Ctrl1", "Ctrl4", "Ctrl5", "Ctrl6", "Ctrl7"]
KEEP_PD = ["PD1", "PD2", "PD6", "PD7", "PD8"]
SAMPLE_RENAME = {sample: sample for sample in KEEP_CTRL + KEEP_PD}
CTRL = KEEP_CTRL
PD = KEEP_PD
SAMPLES = CTRL + PD

STAT_COLS = ["MEAN Ctrl", "MEAN PD", "VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange", "Diff", "Regulation"]
SEARCH_COLS = ["MS2 name", "MS1 name", "Chinese.Name", "HMDB", "KEGG COMPOUND ID", "Formula", "CAS", "INCHIKEY"]

TARGETS = [
    {
        "target": "p-Cresol / 对甲酚 / 4-Cresol",
        "aliases": ["p-cresol", "4-cresol", "para-cresol", "4-methylphenol", "p cresol", "对甲酚"],
        "exclude": ["sulfate", "glucuronide", "acetate", "dinitro", "complexone", "phthalein"],
    },
    {
        "target": "p-Cresol sulfate",
        "aliases": ["p-cresol sulfate", "4-cresol sulfate", "p-cresyl sulfate", "4-cresyl sulfate", "cresol sulfate"],
        "exclude": [],
    },
    {
        "target": "p-Cresol glucuronide",
        "aliases": ["p-cresol glucuronide", "4-cresol glucuronide", "p-cresyl glucuronide", "4-cresyl glucuronide", "cresol glucuronide"],
        "exclude": [],
    },
]


def configure_base_globals() -> None:
    base.RAW_CTRL = RAW_CTRL
    base.RAW_PD = RAW_PD
    base.DROP_RAW = DROP_RAW
    base.KEEP_CTRL = KEEP_CTRL
    base.KEEP_PD = KEEP_PD
    base.SAMPLE_RENAME = SAMPLE_RENAME
    base.CTRL = CTRL
    base.PD = PD
    base.SAMPLES = SAMPLES


def md_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "无。"
    out = frame.copy()
    for col in out.columns:
        if np.issubdtype(out[col].dtype, np.number):
            out[col] = out[col].map(lambda x: f"{x:.4g}" if pd.notna(x) else "")
    out = out.fillna("").astype(str)
    headers = list(out.columns)
    rows = out.values.tolist()
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        safe = [cell.replace("|", "/") for cell in row]
        lines.append("| " + " | ".join(safe) + " |")
    return "\n".join(lines)


def kegg_enrichment(df: pd.DataFrame, source: Path, outdir: Path) -> pd.DataFrame:
    base.ensure_dir(outdir)
    path_file = source / "09.KEGG Analysis" / "KEGG Pathway.xlsx"
    if not path_file.exists():
        return pd.DataFrame()
    path_df = pd.read_excel(path_file, sheet_name="Ctrl vs PD")
    diff_df = df[df["Diff"]].copy()

    all_ids = set()
    for value in df["KEGG COMPOUND ID"]:
        all_ids.update(base.split_kegg_ids(value))
    diff_ids = set()
    up_ids = set()
    down_ids = set()
    for _, row in diff_df.iterrows():
        ids = base.split_kegg_ids(row["KEGG COMPOUND ID"])
        diff_ids.update(ids)
        if row["Log_Foldchange"] > 0:
            up_ids.update(ids)
        elif row["Log_Foldchange"] < 0:
            down_ids.update(ids)
    if not all_ids or not diff_ids:
        return pd.DataFrame()

    m = len(all_ids)
    n = len(diff_ids)
    rows = []
    for _, row in path_df.iterrows():
        pathway_all = set(base.split_kegg_ids(row.get("Compounds (all)", "")))
        if not pathway_all:
            pathway_all = set(base.split_kegg_ids(row.get("Compounds.(all)", "")))
        overlap = pathway_all & diff_ids
        if not overlap:
            continue
        x = len(overlap)
        k = len(pathway_all & all_ids)
        total = row.get("Total", np.nan)
        p_value = hypergeom.sf(x - 1, m, k, n) if k > 0 else np.nan
        rows.append(
            {
                "Pathway": row["Pathway"],
                "Description": str(row["Description"]).replace(" - Mus musculus (house mouse)", ""),
                "#.compounds.(all)": k,
                "Compounds.(all)": ";".join(sorted(pathway_all & all_ids)),
                "#.compounds.(dem)": x,
                "Compounds.(dem)": ";".join(sorted(overlap)),
                "Total": total,
                "Percent": (k / total * 100) if pd.notna(total) and total else np.nan,
                "Rich_factor": x / k if k else np.nan,
                "p_value": p_value,
                "up_nums": len(overlap & up_ids),
                "down_nums": len(overlap & down_ids),
                "DA_score": (len(overlap & up_ids) - len(overlap & down_ids)) / x if x else np.nan,
            }
        )

    enrich = pd.DataFrame(rows)
    if enrich.empty:
        return enrich
    enrich["q_value"] = base.bh_fdr(enrich["p_value"].to_numpy())
    enrich = enrich.sort_values(["p_value", "#.compounds.(dem)"], ascending=[True, False])
    with pd.ExcelWriter(outdir / "KEGG Enrichment data matrix.xlsx") as writer:
        enrich.to_excel(writer, index=False, sheet_name="Ctrl vs PD")

    top = enrich.head(20).copy()
    if not top.empty:
        top = top.sort_values("Rich_factor")
        fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
        sizes = top["#.compounds.(dem)"] * 30
        scatter = ax.scatter(top["Rich_factor"], top["Description"], s=sizes, c=-np.log10(top["p_value"]), cmap="viridis")
        ax.set_xlabel("Rich factor")
        ax.set_ylabel("")
        ax.set_title("KEGG enrichment bubble: Ctrl vs PD")
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label("-log10(P-value)")
        fig.tight_layout()
        fig.savefig(outdir / "KEGG Enrichment bubble.png")
        fig.savefig(outdir / "KEGG Enrichment bubble.pdf")
        plt.close(fig)
    return enrich


def normalize_text(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    text = text.lower().replace("–", "-").replace("—", "-").replace("_", "-")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def row_text(row: pd.Series) -> str:
    return " | ".join(normalize_text(row.get(col, "")) for col in SEARCH_COLS if col in row.index)


def target_mask(df: pd.DataFrame, spec: dict[str, object]) -> pd.Series:
    text = df.apply(row_text, axis=1)
    mask = pd.Series(False, index=df.index)
    for alias in spec["aliases"]:
        mask = mask | text.str.contains(re.escape(normalize_text(alias)), regex=True, na=False)
    for exclude in spec.get("exclude", []):
        mask = mask & ~text.str.contains(re.escape(normalize_text(exclude)), regex=True, na=False)
    return mask


def all_cresol_like_mask(df: pd.DataFrame) -> pd.Series:
    text = df.apply(row_text, axis=1)
    return text.str.contains(r"cresol|cresyl|对甲酚", case=False, regex=True, na=False)


def make_long_table(df: pd.DataFrame) -> pd.DataFrame:
    id_cols = [col for col in ["Target"] + base.META_COLS + STAT_COLS if col in df.columns]
    if df.empty:
        return pd.DataFrame()
    long_df = df.melt(id_vars=id_cols, value_vars=SAMPLES, var_name="Sample", value_name="Intensity")
    long_df["Group"] = long_df["Sample"].map(lambda sample: "Ctrl" if sample in CTRL else "PD")
    keep_cols = ["Target", "ID", "MS2 name", "MS1 name", "HMDB", "Sample", "Group", "Intensity", "MEAN Ctrl", "MEAN PD", "Fold_Change", "Log_Foldchange", "P-Value", "Q-Value", "VIP", "Regulation"]
    return long_df[[col for col in keep_cols if col in long_df.columns]]


def plot_target_boxes(long_df: pd.DataFrame, outdir: Path) -> None:
    if long_df.empty:
        return
    sns.set_theme(style="whitegrid", font="DejaVu Sans")
    for target, target_df in long_df.groupby("Target", sort=False):
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", target).strip("_")
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        sns.boxplot(data=target_df, x="Group", y="Intensity", order=["Ctrl", "PD"], palette=["#4C78A8", "#E45756"], width=0.5, showfliers=False, ax=ax)
        sns.stripplot(data=target_df, x="Group", y="Intensity", order=["Ctrl", "PD"], color="#222222", size=4, jitter=0.14, ax=ax)
        row = target_df.iloc[0]
        title_name = str(row["MS2 name"]) if pd.notna(row.get("MS2 name")) and str(row.get("MS2 name")).strip() else target
        ax.set_title(title_name, fontsize=11, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Intensity")
        ax.text(
            0.5,
            -0.22,
            "FC={0:.3g}, P={1:.3g}, Q={2:.3g}".format(row["Fold_Change"], row["P-Value"], row["Q-Value"]),
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=9,
        )
        for suffix in [".png", ".pdf"]:
            fig.savefig(outdir / f"{safe}_boxplot{suffix}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def export_cresol_targets(stat_table: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    target_dir = outdir / "p_cresol_related_metabolites"
    base.ensure_dir(target_dir)
    ordered_cols = [col for col in base.META_COLS + SAMPLES + STAT_COLS if col in stat_table.columns]
    target_rows = []
    summary_rows = []
    for spec in TARGETS:
        matched = stat_table.loc[target_mask(stat_table, spec), ordered_cols].copy()
        matched.insert(0, "Target", spec["target"])
        target_rows.append(matched)
        summary_rows.append(
            {
                "Target": spec["target"],
                "Found": len(matched) > 0,
                "Matched_rows": len(matched),
                "Matched_IDs": "; ".join(matched["ID"].astype(str).tolist()) if len(matched) else "",
                "Matched_MS2_names": "; ".join(matched["MS2 name"].fillna("").astype(str).tolist()) if len(matched) else "",
                "Note": "" if len(matched) else "No exact target match after removing Ctrl2, Ctrl3, Ctrl8, PD3, PD4, PD5.",
            }
        )

    exact_targets = pd.concat(target_rows, ignore_index=True) if target_rows else pd.DataFrame()
    cresol_like = stat_table.loc[all_cresol_like_mask(stat_table), ordered_cols].copy()
    cresol_like.insert(0, "Target", "all cresol-like matches")
    missing = pd.DataFrame([row for row in summary_rows if not row["Found"]])
    summary = pd.DataFrame(summary_rows)
    long_targets = make_long_table(exact_targets)
    long_cresol_like = make_long_table(cresol_like)

    workbook = target_dir / "midbrain_p_cresol_related_metabolites_after_sample_filter.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="target_summary", index=False)
        exact_targets.to_excel(writer, sheet_name="exact_target_matches", index=False)
        long_targets.to_excel(writer, sheet_name="exact_target_long", index=False)
        cresol_like.to_excel(writer, sheet_name="all_cresol_like_matches", index=False)
        long_cresol_like.to_excel(writer, sheet_name="all_cresol_like_long", index=False)
        missing.to_excel(writer, sheet_name="missing_targets", index=False)

    exact_targets.to_csv(target_dir / "exact_target_matches.csv", index=False, encoding="utf-8-sig")
    long_targets.to_csv(target_dir / "exact_target_sample_values_long.csv", index=False, encoding="utf-8-sig")
    cresol_like.to_csv(target_dir / "all_cresol_like_matches.csv", index=False, encoding="utf-8-sig")
    long_cresol_like.to_csv(target_dir / "all_cresol_like_sample_values_long.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(target_dir / "target_summary.csv", index=False, encoding="utf-8-sig")
    plot_target_boxes(long_targets, target_dir)
    (target_dir / "README.txt").write_text(
        "\n".join(
            [
                "Midbrain p-cresol related metabolites after sample filtering",
                "",
                "Removed samples: Ctrl2, Ctrl3, Ctrl8, PD3, PD4, PD5.",
                "Kept samples: " + ", ".join(SAMPLES),
                "",
                "Exact target results:",
                summary.to_string(index=False),
                "",
                "Broad cresol-like sheet includes all names containing cresol/cresyl/对甲酚.",
            ]
        ),
        encoding="utf-8",
    )
    return summary


def write_report(df: pd.DataFrame, outdir: Path, source: Path, enrich: pd.DataFrame, model_df: pd.DataFrame, target_summary: pd.DataFrame) -> None:
    diff = df[df["Diff"]]
    up = diff[diff["Log_Foldchange"] > 0]
    down = diff[diff["Log_Foldchange"] < 0]
    named_all = df["MS2 name"].dropna().astype(str)
    named_diff = diff["MS2 name"].dropna().astype(str)
    top_up = up.sort_values("P-Value").head(10)[["Feature_Label", "P-Value", "Q-Value", "Fold_Change", "VIP"]]
    top_down = down.sort_values("P-Value").head(10)[["Feature_Label", "P-Value", "Q-Value", "Fold_Change", "VIP"]]
    top_path = enrich.head(10)[["Pathway", "Description", "#.compounds.(dem)", "Rich_factor", "p_value", "q_value"]] if not enrich.empty else pd.DataFrame()

    text = f"""# 中脑非靶代谢组 Ctrl vs PD 重分析

## 分析设定

- 参考原始目录：`{source}`
- 新输出目录：`{outdir}`
- 删除样本：Ctrl2、Ctrl3、Ctrl8、PD3、PD4、PD5
- 保留样本：{", ".join(SAMPLES)}
- 差异筛选口径：沿用 Biotree 常用统计口径，`VIP > 1` 且 `P-Value < 0.05`；P 值为 Bartlett 方差齐性检验后切换等方差 t 检验或 Welch t 检验；Q 值为 BH-FDR。
- 说明：本地复跑使用 PLS-DA 计算 VIP，不能完全等同于原 Biotree OPLS-DA 专有流程，但样本过滤、处理后矩阵和差异统计均从原报告表重新计算。

## 结果概览

| 指标 | 数量 |
|---|---:|
| 总 feature 数 | {len(df)} |
| 注释代谢物名称数 | {named_all.nunique()} |
| 差异 feature 数 | {len(diff)} |
| PD 上调 feature 数 | {len(up)} |
| PD 下调 feature 数 | {len(down)} |
| 差异注释代谢物名称数 | {named_diff.nunique()} |

## 模型信息

{md_table(model_df)}

## 对甲酚相关目标导出

{md_table(target_summary)}

目标代谢物单独保存目录：`p_cresol_related_metabolites/`

## PD 上调 Top 10

{md_table(top_up)}

## PD 下调 Top 10

{md_table(top_down)}

## KEGG 富集 Top 10

{md_table(top_path)}

## 主要输出

- 全量统计表：`02.Statistical Analysis/Statistical Analysis Results.xlsx`
- 差异代谢物表：`02.Statistical Analysis/Differentially Expressed Metabolites.xlsx`
- PCA/PLS-DA/火山图：`02.Statistical Analysis/Ctrl vs PD/`
- 差异热图：`03.Hierarchical Clustering Analysis/Ctrl vs PD/heatmap.png`
- 差异箱线图：`04.Boxplot Analysis/Ctrl vs PD/boxplot_up.png`、`boxplot_down.png`
- KEGG 富集：`10.Enrichment Analysis/Ctrl vs PD/KEGG Enrichment bubble.png`
- ROC：`14.ROC Curve/Ctrl vs PD/`
"""
    (outdir / "Ctrl_vs_PD_reanalysis_summary.md").write_text(text, encoding="utf-8")


def report_only(outdir: Path, source: Path) -> None:
    stat_file = outdir / "02.Statistical Analysis" / "Statistical Analysis Results.xlsx"
    model_file = outdir / "02.Statistical Analysis" / "Models Information.xlsx"
    enrich_file = outdir / "10.Enrichment Analysis" / "Ctrl vs PD" / "KEGG Enrichment data matrix.xlsx"
    df = pd.read_excel(stat_file, sheet_name="Ctrl vs PD")
    df["Feature_Label"] = df.apply(base.feature_label, axis=1).map(lambda x: base.clean_name(x, "feature"))
    enrich = pd.read_excel(enrich_file, sheet_name="Ctrl vs PD") if enrich_file.exists() else pd.DataFrame()
    model_df = pd.read_excel(model_file, sheet_name="Models Information") if model_file.exists() else pd.DataFrame()
    target_summary = export_cresol_targets(df, outdir)
    write_report(df, outdir, source, enrich, model_df, target_summary)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    configure_base_globals()
    source = args.source
    outdir = args.out
    base.ensure_dir(outdir)

    if args.report_only:
        report_only(outdir, source)
        print("Report written:", outdir / "Ctrl_vs_PD_reanalysis_summary.md")
        return

    print("Reading processed Biotree matrix...")
    stat_file = source / "02.Statistical Analysis" / "Statistical Analysis Results.xlsx"
    df = pd.read_excel(stat_file, sheet_name="Ctrl vs PD")
    for col in base.META_COLS:
        if col not in df.columns:
            df[col] = np.nan
    df = df[base.META_COLS + SAMPLES + ["VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange"]]

    print("Recomputing PLS-DA VIP...")
    stat_out = outdir / "02.Statistical Analysis"
    comp_out = stat_out / "Ctrl vs PD"
    base.ensure_dir(comp_out)
    pca_scores = base.plot_pca(df, comp_out)
    pls_scores, vip_df = base.plot_plsda(df, comp_out)
    df = df.drop(columns=["VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange"], errors="ignore")
    df = df.merge(vip_df, on="ID", how="left")

    print("Recomputing univariate statistics...")
    ctrl_arr = df[CTRL].astype(float).to_numpy()
    pd_arr = df[PD].astype(float).to_numpy()
    mean_ctrl = np.nanmean(ctrl_arr, axis=1)
    mean_pd = np.nanmean(pd_arr, axis=1)
    pvals = np.array([base.switched_ttest(ctrl_arr[i], pd_arr[i]) for i in range(len(df))], dtype=float)
    qvals = base.bh_fdr(pvals)
    fc = np.divide(mean_pd, mean_ctrl, out=np.full_like(mean_pd, np.inf), where=mean_ctrl != 0)
    logfc = np.log2(fc, where=fc > 0, out=np.full_like(fc, np.nan))
    df["MEAN Ctrl"] = mean_ctrl
    df["MEAN PD"] = mean_pd
    df["P-Value"] = pvals
    df["Q-Value"] = qvals
    df["Fold_Change"] = fc
    df["Log_Foldchange"] = logfc
    df["Diff"] = (df["VIP"] > 1) & (df["P-Value"] < 0.05)
    df["Regulation"] = np.where(df["Diff"] & (df["Log_Foldchange"] > 0), "up", np.where(df["Diff"] & (df["Log_Foldchange"] < 0), "down", "no_diff"))
    df["Feature_Label"] = df.apply(base.feature_label, axis=1).map(lambda x: base.clean_name(x, "feature"))

    stat_cols = base.META_COLS + CTRL + ["MEAN Ctrl"] + PD + ["MEAN PD", "VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange", "Diff", "Regulation"]
    stat_table = df[stat_cols]
    mean_table = df[base.META_COLS[:9] + CTRL + ["MEAN Ctrl"] + PD + ["MEAN PD"]]
    diff_table = stat_table[stat_table["Diff"]].sort_values("P-Value")

    print("Writing Excel tables...")
    with pd.ExcelWriter(stat_out / "Statistical Analysis Results.xlsx") as writer:
        stat_table.to_excel(writer, index=False, sheet_name="Ctrl vs PD")
    with pd.ExcelWriter(stat_out / "Mean.xlsx") as writer:
        mean_table.to_excel(writer, index=False, sheet_name="Mean")
    feature_summary = pd.DataFrame(
        [
            {
                "Group": "Ctrl vs PD",
                "Feature_all": len(stat_table),
                "Feature_diff": len(diff_table),
                "Feature_diff_up": int((diff_table["Log_Foldchange"] > 0).sum()),
                "Feature_diff_down": int((diff_table["Log_Foldchange"] < 0).sum()),
            }
        ]
    )
    cpd_all = df["MS2 name"].dropna().astype(str).nunique()
    cpd_summary = pd.DataFrame(
        [
            {
                "Group": "Ctrl vs PD",
                "Cpd_all": cpd_all,
                "Cpd_diff": diff_table["MS2 name"].dropna().astype(str).nunique(),
                "Cpd_diff_up": diff_table[diff_table["Log_Foldchange"] > 0]["MS2 name"].dropna().astype(str).nunique(),
                "Cpd_diff_down": diff_table[diff_table["Log_Foldchange"] < 0]["MS2 name"].dropna().astype(str).nunique(),
            }
        ]
    )
    with pd.ExcelWriter(stat_out / "Differentially Expressed Metabolites.xlsx") as writer:
        feature_summary.to_excel(writer, index=False, sheet_name="feature summary")
        cpd_summary.to_excel(writer, index=False, sheet_name="metabolite summary")
        diff_table.to_excel(writer, index=False, sheet_name="Ctrl vs PD")

    pca_r2x = np.nan
    if not pca_scores.empty:
        x, _, _ = base.prepare_matrix(df)
        pca = PCA(n_components=min(3, x.shape[0] - 1, x.shape[1]), random_state=0).fit(x)
        pca_r2x = float(pca.explained_variance_ratio_.sum())
    pls_r2y = float(pls_scores["R2Y"].iloc[0]) if not pls_scores.empty else np.nan
    pls_q2 = float(pls_scores["Q2_LOO"].iloc[0]) if not pls_scores.empty else np.nan
    model_df = pd.DataFrame(
        [
            {"Model": "Model 1", "Type": "PCA", "A": 3, "N": len(SAMPLES), "R2X(cum)": pca_r2x, "R2Y(cum)": np.nan, "Q2(cum)": np.nan, "Title": "Ctrl vs PD"},
            {"Model": "Model 2", "Type": "PLS-DA", "A": 2, "N": len(SAMPLES), "R2X(cum)": np.nan, "R2Y(cum)": pls_r2y, "Q2(cum)": pls_q2, "Title": "Ctrl vs PD"},
        ]
    )
    with pd.ExcelWriter(stat_out / "Models Information.xlsx") as writer:
        model_df.to_excel(writer, index=False, sheet_name="Models Information")

    sample_summary = pd.DataFrame(
        [
            {"Original sample": raw, "New sample": SAMPLE_RENAME.get(raw, ""), "Group": "Ctrl" if raw in KEEP_CTRL else "PD" if raw in KEEP_PD else "Removed"}
            for raw in RAW_CTRL + RAW_PD
        ]
    )
    with pd.ExcelWriter(outdir / "sample_filtering_summary.xlsx") as writer:
        sample_summary.to_excel(writer, index=False, sheet_name="sample filtering")

    print("Exporting p-cresol related metabolites...")
    target_summary = export_cresol_targets(stat_table, outdir)

    print("Drawing plots...")
    base.plot_volcano(df, comp_out)
    heatmap_out = outdir / "03.Hierarchical Clustering Analysis" / "Ctrl vs PD"
    base.ensure_dir(heatmap_out)
    base.plot_heatmap(df, heatmap_out)
    box_out = outdir / "04.Boxplot Analysis" / "Ctrl vs PD"
    base.plot_boxplots(df, box_out)
    corr_out = outdir / "07.Correlation Analysis" / "Ctrl vs PD"
    base.plot_correlation(df, corr_out)
    roc_out = outdir / "14.ROC Curve" / "Ctrl vs PD"
    base.roc_analysis(df, roc_out)

    print("Running KEGG enrichment...")
    enrich_out = outdir / "10.Enrichment Analysis" / "Ctrl vs PD"
    enrich = kegg_enrichment(df, source, enrich_out)

    print("Writing report...")
    write_report(df, outdir, source, enrich, model_df, target_summary)
    print("Done:", outdir)


if __name__ == "__main__":
    main()
