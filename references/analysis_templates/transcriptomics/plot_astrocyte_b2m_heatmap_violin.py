from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


matplotlib.use("Agg")

OUT = Path("p_cresol_snRNA_analysis")
FIG = OUT / "figures"
TAB = OUT / "tables"
H5 = OUT / "h5ad"


def save_df(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def p_to_stars(p: float) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def get_gene_expr(adata: ad.AnnData, gene: str) -> np.ndarray:
    source = adata.raw.to_adata() if adata.raw is not None else adata
    if gene not in source.var_names:
        raise ValueError(f"{gene} not found in AnnData.var_names")
    x = source[:, [gene]].X
    return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()


def summarize_b2m(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subtype, sub in df.groupby("astrocyte_refined_subtype", observed=False):
        ctrl = sub.loc[sub["group"] == "Ctrl", "B2m"]
        pdv = sub.loc[sub["group"] == "PD", "B2m"]
        p = stats.mannwhitneyu(pdv, ctrl, alternative="two-sided").pvalue if len(ctrl) and len(pdv) else np.nan
        p_greater = stats.mannwhitneyu(pdv, ctrl, alternative="greater").pvalue if len(ctrl) and len(pdv) else np.nan
        if len(ctrl) and len(pdv):
            rank_stat, rank_p = stats.ranksums(pdv, ctrl)
            rank_p_greater = rank_p / 2 if rank_stat > 0 else 1 - rank_p / 2
        else:
            rank_stat, rank_p, rank_p_greater = np.nan, np.nan, np.nan
        rows.append(
            {
                "astrocyte_refined_subtype": subtype,
                "n_Ctrl": int(len(ctrl)),
                "n_PD": int(len(pdv)),
                "Ctrl_mean_B2m": float(ctrl.mean()) if len(ctrl) else np.nan,
                "PD_mean_B2m": float(pdv.mean()) if len(pdv) else np.nan,
                "PD_minus_Ctrl": float(pdv.mean() - ctrl.mean()) if len(ctrl) and len(pdv) else np.nan,
                "Ctrl_pct_B2m_pos": float((ctrl > 0).mean() * 100) if len(ctrl) else np.nan,
                "PD_pct_B2m_pos": float((pdv > 0).mean() * 100) if len(pdv) else np.nan,
                "cell_level_mwu_p_descriptive": float(p) if np.isfinite(p) else np.nan,
                "cell_level_mwu_p_PD_greater_exploratory": float(p_greater) if np.isfinite(p_greater) else np.nan,
                "cell_level_ranksum_stat": float(rank_stat) if np.isfinite(rank_stat) else np.nan,
                "cell_level_ranksum_p_two_sided": float(rank_p) if np.isfinite(rank_p) else np.nan,
                "cell_level_ranksum_p_PD_greater": float(rank_p_greater) if np.isfinite(rank_p_greater) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    valid = out["cell_level_mwu_p_descriptive"].notna()
    out["cell_level_mwu_fdr_descriptive"] = np.nan
    if valid.any():
        out.loc[valid, "cell_level_mwu_fdr_descriptive"] = multipletests(
            out.loc[valid, "cell_level_mwu_p_descriptive"], method="fdr_bh"
        )[1]
    out["stars"] = out["cell_level_mwu_p_descriptive"].map(p_to_stars)
    out["stars_PD_greater_exploratory"] = out["cell_level_mwu_p_PD_greater_exploratory"].map(p_to_stars)
    out["stars_ranksum_two_sided"] = out["cell_level_ranksum_p_two_sided"].map(p_to_stars)
    out["stars_ranksum_PD_greater"] = out["cell_level_ranksum_p_PD_greater"].map(p_to_stars)
    return out


def plot_heatmap(summary: pd.DataFrame, order: list[str]) -> None:
    mean_mat = summary.set_index("astrocyte_refined_subtype")[["Ctrl_mean_B2m", "PD_mean_B2m"]].loc[order]
    mean_mat.columns = ["Ctrl", "PD"]
    delta = summary.set_index("astrocyte_refined_subtype").loc[order, "PD_minus_Ctrl"]
    stars = summary.set_index("astrocyte_refined_subtype").loc[order, "stars"]

    annot = mean_mat.copy().astype(str)
    for subtype in order:
        annot.loc[subtype, "Ctrl"] = f"{mean_mat.loc[subtype, 'Ctrl']:.2f}"
        annot.loc[subtype, "PD"] = f"{mean_mat.loc[subtype, 'PD']:.2f}\n{stars.loc[subtype]}"

    fig, ax = plt.subplots(figsize=(7.0, max(4.8, 0.46 * len(order))))
    sns.heatmap(
        mean_mat,
        cmap="Reds",
        annot=annot,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Mean B2m expression"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Astrocyte subtype B2m expression")
    fig.subplots_adjust(left=0.52, right=0.92, top=0.92, bottom=0.06)
    fig.savefig(FIG / "astrocyte_subtype_B2m_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    delta_df = delta.to_frame("PD_minus_Ctrl")
    fig, ax = plt.subplots(figsize=(6.0, max(4.8, 0.46 * len(order))))
    sns.heatmap(
        delta_df,
        cmap="vlag",
        center=0,
        annot=delta_df.map(lambda x: f"{x:+.2f}"),
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "PD - Ctrl"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("B2m change")
    fig.subplots_adjust(left=0.58, right=0.88, top=0.92, bottom=0.06)
    fig.savefig(FIG / "astrocyte_subtype_B2m_PD_minus_Ctrl_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_violin(
    df: pd.DataFrame,
    summary: pd.DataFrame,
    order: list[str],
    filename: str,
    star_col: str,
    title: str,
    footnote: str,
) -> None:
    fig_w = max(12, 0.82 * len(order))
    fig, ax = plt.subplots(figsize=(fig_w, 5.8))
    sns.violinplot(
        data=df,
        x="astrocyte_refined_subtype",
        y="B2m",
        hue="group",
        order=order,
        hue_order=["Ctrl", "PD"],
        split=False,
        inner=None,
        cut=0,
        linewidth=0.6,
        palette={"Ctrl": "#4f79a7", "PD": "#c45a5a"},
        ax=ax,
    )
    sns.stripplot(
        data=df.sample(min(len(df), 1800), random_state=1337),
        x="astrocyte_refined_subtype",
        y="B2m",
        hue="group",
        order=order,
        hue_order=["Ctrl", "PD"],
        dodge=True,
        size=1.2,
        alpha=0.18,
        linewidth=0,
        palette={"Ctrl": "#1f4e79", "PD": "#933333"},
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    keep = []
    seen = set()
    for h, l in zip(handles, labels):
        if l in {"Ctrl", "PD"} and l not in seen:
            keep.append((h, l))
            seen.add(l)
    ax.legend([h for h, _ in keep], [l for _, l in keep], title="Group", frameon=False, loc="upper right")

    ymax = float(df["B2m"].quantile(0.995))
    if ymax <= 0:
        ymax = float(df["B2m"].max()) + 0.2
    y_pad = max(0.12, ymax * 0.08)
    stats_by = summary.set_index("astrocyte_refined_subtype")
    for i, subtype in enumerate(order):
        if subtype not in stats_by.index:
            continue
        sub = df[df["astrocyte_refined_subtype"] == subtype]
        local_max = float(sub["B2m"].quantile(0.995)) if len(sub) else ymax
        y = local_max + y_pad
        ax.plot([i - 0.22, i - 0.22, i + 0.22, i + 0.22], [y, y + y_pad * 0.25, y + y_pad * 0.25, y], color="black", lw=0.8)
        ax.text(i, y + y_pad * 0.32, stats_by.loc[subtype, star_col], ha="center", va="bottom", fontsize=9)

    ax.set_ylim(0, max(ymax + y_pad * 3.2, ax.get_ylim()[1]))
    ax.set_xlabel("")
    ax.set_ylabel("B2m expression (log-normalized)")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        0.01,
        footnote,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIG / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    adata = ad.read_h5ad(H5 / "astrocyte_subclustered_refined.h5ad")
    df = adata.obs[["group", "astrocyte_refined_subtype"]].copy()
    df["B2m"] = get_gene_expr(adata, "B2m")
    df["astrocyte_refined_subtype"] = df["astrocyte_refined_subtype"].astype(str)
    summary = summarize_b2m(df)
    summary = summary.sort_values(["PD_minus_Ctrl", "PD_mean_B2m"], ascending=False)
    order = summary["astrocyte_refined_subtype"].tolist()
    save_df(summary, TAB / "astrocyte_subtype_B2m_violin_heatmap_stats.csv")
    plot_heatmap(summary, order)
    plot_violin(
        df,
        summary,
        order,
        "astrocyte_subtype_B2m_violin_with_stats.png",
        "stars",
        "B2m expression across astrocyte refined subtypes",
        "Stars: cell-level two-sided Mann-Whitney U test; descriptive only because Ctrl/PD each has one library.",
    )
    plot_violin(
        df,
        summary,
        order,
        "astrocyte_subtype_B2m_violin_PD_greater_exploratory.png",
        "stars_PD_greater_exploratory",
        "B2m expression across astrocyte refined subtypes (PD > Ctrl exploratory)",
        "Stars: exploratory one-sided Mann-Whitney U test for PD > Ctrl; unadjusted and descriptive only.",
    )
    plot_violin(
        df,
        summary,
        order,
        "astrocyte_subtype_B2m_violin_ranksum_PD_greater.png",
        "stars_ranksum_PD_greater",
        "B2m expression across astrocyte refined subtypes (rank-sum PD > Ctrl)",
        "Stars: exploratory one-sided Wilcoxon rank-sum test for PD > Ctrl; unadjusted and descriptive only.",
    )
    (OUT / "ASTRO_B2M_热图小提琴图说明.md").write_text(
        "\n".join(
            [
                "# 星胶亚群 B2m 表达图",
                "",
                "新增输出：",
                "- `figures/astrocyte_subtype_B2m_heatmap.png`",
                "- `figures/astrocyte_subtype_B2m_PD_minus_Ctrl_heatmap.png`",
                "- `figures/astrocyte_subtype_B2m_violin_with_stats.png`",
                "- `figures/astrocyte_subtype_B2m_violin_PD_greater_exploratory.png`",
                "- `figures/astrocyte_subtype_B2m_violin_ranksum_PD_greater.png`",
                "- `tables/astrocyte_subtype_B2m_violin_heatmap_stats.csv`",
                "",
                "小提琴图星号为细胞级 Mann-Whitney U test：ns, *, **, ***, ****。",
                "新增 `PD_greater_exploratory` 图使用单侧 Mann-Whitney U test，方向为 PD > Ctrl，未做多重校正，仅适合有明确先验方向假设时作为探索性补充。",
                "新增 `ranksum_PD_greater` 图使用单侧 Wilcoxon rank-sum test，方向为 PD > Ctrl，未做多重校正。",
                "由于 Ctrl/PD 各一个 library，该检验仅用于描述性标注，不应作为动物重复层面的显著性结论。",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
