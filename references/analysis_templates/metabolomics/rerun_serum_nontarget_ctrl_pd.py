from __future__ import annotations

import argparse
import math
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import hypergeom
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler


META_COLS = [
    "ID",
    "MS2 name",
    "MS2 score",
    "level",
    "MS1 name",
    "MS1 ppm",
    "mz",
    "rt",
    "type",
    "Formula",
    "HMDB",
    "Chinese.Name",
    "CAS",
    "KEGG COMPOUND ID",
    "Super.Class",
    "Superclass.Chinese",
    "Class",
    "Class.Chinese",
    "Sub.Class",
    "Subclass.Chinese",
    "INCHIKEY",
]

RAW_CTRL = [f"SC{i}" for i in range(1, 9)]
RAW_PD = [f"SP{i}" for i in range(1, 9)]
DROP_RAW = ["SC1", "SC6", "SP4", "SP6"]
KEEP_CTRL = ["SC2", "SC3", "SC4", "SC5", "SC7", "SC8"]
KEEP_PD = ["SP1", "SP2", "SP3", "SP5", "SP7", "SP8"]
SAMPLE_RENAME = {
    "SC2": "Ctrl2",
    "SC3": "Ctrl3",
    "SC4": "Ctrl4",
    "SC5": "Ctrl5",
    "SC7": "Ctrl7",
    "SC8": "Ctrl8",
    "SP1": "PD1",
    "SP2": "PD2",
    "SP3": "PD3",
    "SP5": "PD5",
    "SP7": "PD7",
    "SP8": "PD8",
}
CTRL = [SAMPLE_RENAME[c] for c in KEEP_CTRL]
PD = [SAMPLE_RENAME[c] for c in KEEP_PD]
SAMPLES = CTRL + PD


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


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


def switched_ttest(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if np.nanstd(a) == 0 and np.nanstd(b) == 0:
        return 1.0 if np.nanmean(a) == np.nanmean(b) else 0.0
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p_var = stats.bartlett(a, b).pvalue
        equal_var = bool(np.isfinite(p_var) and p_var >= 0.05)
    except Exception:
        equal_var = False
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = stats.ttest_ind(a, b, equal_var=equal_var, nan_policy="omit").pvalue
        if not np.isfinite(p):
            return 1.0
        return float(p)
    except Exception:
        return np.nan


def pls_vip(pls: PLSRegression) -> np.ndarray:
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_
    p = w.shape[0]
    ss = np.diag(t.T @ t @ q.T @ q).reshape(-1)
    denom = ss.sum()
    if denom <= 0 or not np.isfinite(denom):
        return np.ones(p)
    weight = (w / np.linalg.norm(w, axis=0)) ** 2
    vip = np.sqrt(p * (weight @ ss) / denom)
    return np.asarray(vip, dtype=float)


def clean_name(value: object, fallback: str) -> str:
    text = "" if pd.isna(value) else str(value)
    if not text.strip():
        text = fallback
    text = re.sub(r"[\\/:*?\"<>|]+", "_", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text[:120] or fallback


def split_kegg_ids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    ids = re.split(r"[;,\s]+", str(value))
    return [x.strip() for x in ids if re.match(r"^C\d{5}$", x.strip())]


def zscore_rows(values: pd.DataFrame) -> pd.DataFrame:
    arr = values.astype(float).to_numpy()
    mean = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    z = (arr - mean) / sd
    return pd.DataFrame(z, index=values.index, columns=values.columns)


def feature_label(row: pd.Series) -> str:
    name = row.get("MS2 name")
    if pd.isna(name) or str(name).strip() == "":
        return f"feature_{int(row['ID'])}" if pd.notna(row.get("ID")) else "feature"
    return str(name)


def prepare_matrix(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    values = df[SAMPLES].astype(float).T
    values = values.replace([np.inf, -np.inf], np.nan).fillna(0)
    eps = max(float(values[values > 0].min().min()) / 2, 1e-12) if (values > 0).any().any() else 1e-12
    log_values = np.log10(values + eps)
    keep_var = log_values.var(axis=0) > 0
    log_values = log_values.loc[:, keep_var]
    scaler = StandardScaler()
    x = scaler.fit_transform(log_values)
    y = np.array([0] * len(CTRL) + [1] * len(PD), dtype=float)
    return x, y, log_values


def plot_pca(df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    x, y, _ = prepare_matrix(df)
    pca = PCA(n_components=min(3, x.shape[0] - 1, x.shape[1]), random_state=0)
    scores = pca.fit_transform(x)
    score_df = pd.DataFrame(
        {
            "Sample": SAMPLES,
            "Group": ["Ctrl"] * len(CTRL) + ["PD"] * len(PD),
            "PC1": scores[:, 0],
            "PC2": scores[:, 1] if scores.shape[1] > 1 else 0,
            "PC3": scores[:, 2] if scores.shape[1] > 2 else 0,
        }
    )
    ensure_dir(outdir)
    fig, ax = plt.subplots(figsize=(6.5, 5.2), dpi=160)
    sns.scatterplot(data=score_df, x="PC1", y="PC2", hue="Group", style="Group", s=85, ax=ax)
    for _, r in score_df.iterrows():
        ax.text(r["PC1"], r["PC2"], r["Sample"], fontsize=8, ha="left", va="bottom")
    ax.axhline(0, color="#cccccc", lw=0.8)
    ax.axvline(0, color="#cccccc", lw=0.8)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    pc2 = pca.explained_variance_ratio_[1] * 100 if len(pca.explained_variance_ratio_) > 1 else 0
    ax.set_ylabel(f"PC2 ({pc2:.1f}%)")
    ax.set_title("PCA score plot: Ctrl vs PD")
    fig.tight_layout()
    fig.savefig(outdir / "PCA score plot.png")
    fig.savefig(outdir / "PCA score plot.pdf")
    plt.close(fig)
    return score_df


def plot_plsda(df: pd.DataFrame, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    x, y, log_values = prepare_matrix(df)
    n_components = min(2, x.shape[0] - 2, x.shape[1])
    if n_components < 1:
        vip = pd.DataFrame({"ID": df["ID"], "VIP": 1.0})
        return pd.DataFrame(), vip
    pls = PLSRegression(n_components=n_components, scale=False)
    pls.fit(x, y)
    scores = pls.x_scores_
    y_pred = pls.predict(x).ravel()
    r2y = 1 - np.sum((y - y_pred) ** 2) / np.sum((y - y.mean()) ** 2)
    try:
        loo_pred = cross_val_predict(pls, x, y, cv=LeaveOneOut()).ravel()
        q2 = 1 - np.sum((y - loo_pred) ** 2) / np.sum((y - y.mean()) ** 2)
    except Exception:
        q2 = np.nan
    score_df = pd.DataFrame(
        {
            "Sample": SAMPLES,
            "Group": ["Ctrl"] * len(CTRL) + ["PD"] * len(PD),
            "PLS1": scores[:, 0],
            "PLS2": scores[:, 1] if scores.shape[1] > 1 else 0,
            "Y_pred": y_pred,
            "R2Y": r2y,
            "Q2_LOO": q2,
        }
    )
    fig, ax = plt.subplots(figsize=(6.5, 5.2), dpi=160)
    sns.scatterplot(data=score_df, x="PLS1", y="PLS2", hue="Group", style="Group", s=85, ax=ax)
    for _, r in score_df.iterrows():
        ax.text(r["PLS1"], r["PLS2"], r["Sample"], fontsize=8, ha="left", va="bottom")
    ax.axhline(0, color="#cccccc", lw=0.8)
    ax.axvline(0, color="#cccccc", lw=0.8)
    ax.set_title(f"PLS-DA score plot: Ctrl vs PD (R2Y={r2y:.3f}, Q2={q2:.3f})")
    fig.tight_layout()
    fig.savefig(outdir / "PLS-DA score plot.png")
    fig.savefig(outdir / "PLS-DA score plot.pdf")
    plt.close(fig)

    vip_values = pls_vip(pls)
    full_vip = pd.Series(1.0, index=df.index)
    full_vip.loc[log_values.columns] = vip_values
    vip = pd.DataFrame({"ID": df["ID"], "VIP": full_vip.to_numpy()})
    return score_df, vip


def plot_volcano(df: pd.DataFrame, outdir: Path) -> None:
    plot_df = df.copy()
    plot_df["neg_log10_p"] = -np.log10(plot_df["P-Value"].clip(lower=1e-300))
    plot_df["status"] = "Not significant"
    plot_df.loc[(plot_df["Diff"] == True) & (plot_df["Log_Foldchange"] > 0), "status"] = "Up in PD"
    plot_df.loc[(plot_df["Diff"] == True) & (plot_df["Log_Foldchange"] < 0), "status"] = "Down in PD"
    palette = {"Not significant": "#bdbdbd", "Up in PD": "#d73027", "Down in PD": "#4575b4"}
    fig, ax = plt.subplots(figsize=(7, 5.4), dpi=160)
    sns.scatterplot(
        data=plot_df,
        x="Log_Foldchange",
        y="neg_log10_p",
        hue="status",
        palette=palette,
        s=8,
        linewidth=0,
        alpha=0.75,
        ax=ax,
    )
    ax.axhline(-math.log10(0.05), color="#555555", lw=0.8, ls="--")
    ax.axvline(0, color="#555555", lw=0.8)
    top = pd.concat(
        [
            plot_df[(plot_df["Diff"]) & (plot_df["Log_Foldchange"] > 0)].nsmallest(10, "P-Value"),
            plot_df[(plot_df["Diff"]) & (plot_df["Log_Foldchange"] < 0)].nsmallest(10, "P-Value"),
        ]
    )
    for _, r in top.iterrows():
        ax.text(r["Log_Foldchange"], r["neg_log10_p"], clean_name(r["MS2 name"], f"feature_{r['ID']}"), fontsize=6)
    ax.set_xlabel("log2(Fold change: PD / Ctrl)")
    ax.set_ylabel("-log10(P-value)")
    ax.set_title("Volcano plot: Ctrl vs PD")
    fig.tight_layout()
    fig.savefig(outdir / "volcano plot.png")
    fig.savefig(outdir / "volcano plot.pdf")
    plt.close(fig)


def plot_heatmap(df: pd.DataFrame, outdir: Path) -> None:
    diff = df[df["Diff"]].sort_values("P-Value").head(50).copy()
    if diff.empty:
        return
    mat = diff.set_index("Feature_Label")[SAMPLES]
    z = zscore_rows(np.log10(mat.astype(float) + 1e-12))
    colors = pd.Series(["#4c78a8"] * len(CTRL) + ["#e45756"] * len(PD), index=SAMPLES)
    g = sns.clustermap(
        z,
        cmap="vlag",
        col_colors=colors,
        figsize=(8, max(8, len(z) * 0.16)),
        xticklabels=True,
        yticklabels=True,
        center=0,
    )
    g.fig.suptitle("Top differential metabolites heatmap: Ctrl vs PD", y=1.02)
    g.savefig(outdir / "heatmap.png", dpi=180, bbox_inches="tight")
    g.savefig(outdir / "heatmap.pdf", bbox_inches="tight")
    plt.close(g.fig)


def plot_boxplots(df: pd.DataFrame, outdir: Path) -> None:
    ensure_dir(outdir)
    for direction, title in [("up", "Top up metabolites in PD"), ("down", "Top down metabolites in PD")]:
        if direction == "up":
            sub = df[(df["Diff"]) & (df["Log_Foldchange"] > 0)].sort_values("P-Value").head(10)
        else:
            sub = df[(df["Diff"]) & (df["Log_Foldchange"] < 0)].sort_values("P-Value").head(10)
        if sub.empty:
            continue
        long = sub[["Feature_Label"] + SAMPLES].melt(id_vars="Feature_Label", var_name="Sample", value_name="Abundance")
        long["Group"] = long["Sample"].map({s: "Ctrl" for s in CTRL} | {s: "PD" for s in PD})
        long["log10 abundance"] = np.log10(long["Abundance"].astype(float) + 1e-12)
        fig, ax = plt.subplots(figsize=(10, 5.8), dpi=160)
        sns.boxplot(data=long, x="Feature_Label", y="log10 abundance", hue="Group", ax=ax, fliersize=2)
        sns.stripplot(
            data=long,
            x="Feature_Label",
            y="log10 abundance",
            hue="Group",
            dodge=True,
            ax=ax,
            size=3,
            alpha=0.7,
            legend=False,
        )
        ax.set_title(title)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=55, labelsize=7)
        fig.tight_layout()
        fig.savefig(outdir / f"boxplot_{direction}.png")
        fig.savefig(outdir / f"boxplot_{direction}.pdf")
        plt.close(fig)


def plot_correlation(df: pd.DataFrame, outdir: Path) -> None:
    ensure_dir(outdir)
    sub = df[df["Diff"]].sort_values("P-Value").head(30)
    if len(sub) < 3:
        return
    mat = sub.set_index("Feature_Label")[SAMPLES].T
    corr = np.log10(mat.astype(float) + 1e-12).corr(method="spearman")
    fig, ax = plt.subplots(figsize=(9, 7.5), dpi=160)
    sns.heatmap(corr, cmap="vlag", center=0, ax=ax, xticklabels=False, yticklabels=False)
    ax.set_title("Spearman correlation of top differential metabolites")
    fig.tight_layout()
    fig.savefig(outdir / "Correlation plot.png")
    fig.savefig(outdir / "Correlation plot.pdf")
    plt.close(fig)
    corr.to_excel(outdir / "correlational matrix.xlsx")


def kegg_enrichment(df: pd.DataFrame, source: Path, outdir: Path) -> pd.DataFrame:
    ensure_dir(outdir)
    path_file = source / "09.KEGG Analysis" / "KEGG Pathway.xlsx"
    if not path_file.exists():
        return pd.DataFrame()
    path_df = pd.read_excel(path_file, sheet_name="Ctrl vs ppm30")
    diff_df = df[df["Diff"]].copy()
    all_ids = set()
    for value in df["KEGG COMPOUND ID"]:
        all_ids.update(split_kegg_ids(value))
    diff_ids = set()
    up_ids = set()
    down_ids = set()
    for _, row in diff_df.iterrows():
        ids = split_kegg_ids(row["KEGG COMPOUND ID"])
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
        pathway_all = set(split_kegg_ids(row.get("Compounds (all)", "")))
        if not pathway_all:
            pathway_all = set(split_kegg_ids(row.get("Compounds.(all)", "")))
        overlap = pathway_all & diff_ids
        if not overlap:
            continue
        x = len(overlap)
        k = len(pathway_all & all_ids)
        total = row.get("Total", np.nan)
        p = hypergeom.sf(x - 1, m, k, n) if k > 0 else np.nan
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
                "p_value": p,
                "up_nums": len(overlap & up_ids),
                "down_nums": len(overlap & down_ids),
                "DA_score": (len(overlap & up_ids) - len(overlap & down_ids)) / x if x else np.nan,
            }
        )
    enrich = pd.DataFrame(rows)
    if enrich.empty:
        return enrich
    enrich["q_value"] = bh_fdr(enrich["p_value"].to_numpy())
    enrich = enrich.sort_values(["p_value", "#.compounds.(dem)"], ascending=[True, False])
    with pd.ExcelWriter(outdir / "KEGG Enrichment data matrix.xlsx") as writer:
        enrich.to_excel(writer, index=False, sheet_name="Ctrl vs PD")
    top = enrich.head(20).copy()
    if not top.empty:
        top = top.sort_values("Rich_factor")
        fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
        size = top["#.compounds.(dem)"] * 30
        sc = ax.scatter(top["Rich_factor"], top["Description"], s=size, c=-np.log10(top["p_value"]), cmap="viridis")
        ax.set_xlabel("Rich factor")
        ax.set_ylabel("")
        ax.set_title("KEGG enrichment bubble: Ctrl vs PD")
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label("-log10(P-value)")
        fig.tight_layout()
        fig.savefig(outdir / "KEGG Enrichment bubble.png")
        fig.savefig(outdir / "KEGG Enrichment bubble.pdf")
        plt.close(fig)
    return enrich


def roc_analysis(df: pd.DataFrame, outdir: Path) -> None:
    ensure_dir(outdir)
    y = np.array([0] * len(CTRL) + [1] * len(PD), dtype=int)
    for direction, filename in [("up", "auc-up.csv"), ("down", "auc-down.csv")]:
        if direction == "up":
            sub = df[(df["Diff"]) & (df["Log_Foldchange"] > 0)].sort_values("P-Value").head(10)
        else:
            sub = df[(df["Diff"]) & (df["Log_Foldchange"] < 0)].sort_values("P-Value").head(10)
        rows = []
        for _, row in sub.iterrows():
            scores = row[SAMPLES].astype(float).to_numpy()
            if direction == "down":
                scores = -scores
            fpr, tpr, thresholds = roc_curve(y, scores)
            auc_value = auc(fpr, tpr)
            rows.append({"Feature_Label": row["Feature_Label"], "AUC": auc_value, "P-Value": row["P-Value"], "Log2FC": row["Log_Foldchange"]})
            fig, ax = plt.subplots(figsize=(4.8, 4.5), dpi=160)
            ax.plot(fpr, tpr, color="#d73027", lw=2, label=f"AUC={auc_value:.3f}")
            ax.plot([0, 1], [0, 1], color="#999999", ls="--", lw=1)
            ax.set_xlabel("False positive rate")
            ax.set_ylabel("True positive rate")
            ax.set_title(row["Feature_Label"][:80])
            ax.legend(loc="lower right")
            fig.tight_layout()
            fig.savefig(outdir / f"{clean_name(row['Feature_Label'], str(row['ID']))}.ROC.png")
            plt.close(fig)
        pd.DataFrame(rows).to_csv(outdir / filename, index=False, encoding="utf-8-sig")


def write_report(df: pd.DataFrame, outdir: Path, enrich: pd.DataFrame, model_df: pd.DataFrame) -> None:
    diff = df[df["Diff"]]
    up = diff[diff["Log_Foldchange"] > 0]
    down = diff[diff["Log_Foldchange"] < 0]
    named_all = df["MS2 name"].dropna().astype(str)
    named_diff = diff["MS2 name"].dropna().astype(str)
    top_up = up.sort_values("P-Value").head(10)[["Feature_Label", "P-Value", "Q-Value", "Fold_Change", "VIP"]]
    top_down = down.sort_values("P-Value").head(10)[["Feature_Label", "P-Value", "Q-Value", "Fold_Change", "VIP"]]
    top_path = enrich.head(10)[["Pathway", "Description", "#.compounds.(dem)", "Rich_factor", "p_value", "q_value"]] if not enrich.empty else pd.DataFrame()

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

    text = f"""# 血清非靶代谢组 Ctrl vs PD 重分析

## 分析设定

- 参考原始目录：`{outdir.parent / '血清非靶'}`
- 新输出目录：`{outdir}`
- 删除样本：Ctrl1/SC1、Ctrl6/SC6、30ppm4/SP4、30ppm6/SP6
- 保留样本：Ctrl2、Ctrl3、Ctrl4、Ctrl5、Ctrl7、Ctrl8；PD1、PD2、PD3、PD5、PD7、PD8
- 分组改名：原 `30ppm`/`SP` 组在本次输出中统一命名为 `PD`
- 差异筛选口径：沿用 Biotree 常用统计口径，`VIP > 1` 且 `P-Value < 0.05`；P 值为 Bartlett 方差齐性检验后切换等方差 t 检验或 Welch t 检验；Q 值为 BH-FDR。
- 说明：本地复跑使用 PLS-DA 计算 VIP，不能完全等同于原 Biotree OPLS-DA 专有流程，但样本过滤、归一化矩阵和差异统计均从原报告处理后矩阵重新计算。

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


def report_only(outdir: Path) -> None:
    stat_file = outdir / "02.Statistical Analysis" / "Statistical Analysis Results.xlsx"
    model_file = outdir / "02.Statistical Analysis" / "Models Information.xlsx"
    enrich_file = outdir / "10.Enrichment Analysis" / "Ctrl vs PD" / "KEGG Enrichment data matrix.xlsx"
    df = pd.read_excel(stat_file, sheet_name="Ctrl vs PD")
    df["Feature_Label"] = df.apply(feature_label, axis=1).map(lambda x: clean_name(x, "feature"))
    enrich = pd.read_excel(enrich_file, sheet_name="Ctrl vs PD") if enrich_file.exists() else pd.DataFrame()
    model_df = pd.read_excel(model_file, sheet_name="Models Information") if model_file.exists() else pd.DataFrame()
    write_report(df, outdir, enrich, model_df)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    source = args.source
    outdir = args.out
    ensure_dir(outdir)

    if args.report_only:
        report_only(outdir)
        print("Report written:", outdir / "Ctrl_vs_PD_reanalysis_summary.md")
        return

    print("Reading processed Biotree matrix...")
    stat_file = source / "02.Statistical Analysis" / "Statistical Analysis Results.xlsx"
    df = pd.read_excel(stat_file, sheet_name="Ctrl vs ppm30")
    df = df.rename(columns=SAMPLE_RENAME)
    for col in META_COLS:
        if col not in df.columns:
            df[col] = np.nan
    df = df[META_COLS + SAMPLES + ["VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange"]]

    print("Recomputing PLS-DA VIP...")
    stat_out = outdir / "02.Statistical Analysis"
    comp_out = stat_out / "Ctrl vs PD"
    ensure_dir(comp_out)
    pca_scores = plot_pca(df, comp_out)
    pls_scores, vip_df = plot_plsda(df, comp_out)
    df = df.drop(columns=["VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange"], errors="ignore")
    df = df.merge(vip_df, on="ID", how="left")

    print("Recomputing univariate statistics...")
    ctrl_arr = df[CTRL].astype(float).to_numpy()
    pd_arr = df[PD].astype(float).to_numpy()
    mean_ctrl = np.nanmean(ctrl_arr, axis=1)
    mean_pd = np.nanmean(pd_arr, axis=1)
    pvals = np.array([switched_ttest(ctrl_arr[i], pd_arr[i]) for i in range(len(df))], dtype=float)
    qvals = bh_fdr(pvals)
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
    df["Feature_Label"] = df.apply(feature_label, axis=1).map(lambda x: clean_name(x, "feature"))

    stat_cols = META_COLS + CTRL + ["MEAN Ctrl"] + PD + ["MEAN PD", "VIP", "P-Value", "Q-Value", "Fold_Change", "Log_Foldchange", "Diff", "Regulation"]
    stat_table = df[stat_cols]
    mean_table = df[META_COLS[:9] + CTRL + ["MEAN Ctrl"] + PD + ["MEAN PD"]]
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
    cpd_diff = diff_table["MS2 name"].dropna().astype(str)
    cpd_summary = pd.DataFrame(
        [
            {
                "Group": "Ctrl vs PD",
                "Cpd_all": cpd_all,
                "Cpd_diff": cpd_diff.nunique(),
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
        # R2X summary is approximated from the saved PCA model by recomputing.
        x, _, _ = prepare_matrix(df)
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

    print("Drawing plots...")
    plot_volcano(df, comp_out)
    heatmap_out = outdir / "03.Hierarchical Clustering Analysis" / "Ctrl vs PD"
    ensure_dir(heatmap_out)
    plot_heatmap(df, heatmap_out)
    box_out = outdir / "04.Boxplot Analysis" / "Ctrl vs PD"
    plot_boxplots(df, box_out)
    corr_out = outdir / "07.Correlation Analysis" / "Ctrl vs PD"
    plot_correlation(df, corr_out)
    roc_out = outdir / "14.ROC Curve" / "Ctrl vs PD"
    roc_analysis(df, roc_out)

    print("Running KEGG enrichment...")
    enrich_out = outdir / "10.Enrichment Analysis" / "Ctrl vs PD"
    enrich = kegg_enrichment(df, source, enrich_out)

    print("Writing report...")
    write_report(df, outdir, enrich, model_df)
    print("Done:", outdir)


if __name__ == "__main__":
    main()
