from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

from plot_astrocyte_refined_c6_c12_merged_embeddings import RENAMED_COL, add_renamed_label


matplotlib.use("Agg")
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=True)

OUT = Path("p_cresol_snRNA_analysis")
FIG = OUT / "figures"
TAB = OUT / "tables"
H5 = OUT / "h5ad"

GROUPBY = RENAMED_COL
RANK_KEY = "rank_genes_groups_astrocyte_c6_c12_merged_top_markers"
TOP_N_PER_SUBTYPE = 3

EXCLUDE_EXACT = {"Malat1", "Xist"}
EXCLUDE_PREFIXES = ("mt-", "Rpl", "Rps", "Gm")
EXCLUDE_SUFFIXES = ("Rik",)


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


def is_informative_gene(gene: str) -> bool:
    if gene in EXCLUDE_EXACT:
        return False
    if any(gene.startswith(prefix) for prefix in EXCLUDE_PREFIXES):
        return False
    if any(gene.endswith(suffix) for suffix in EXCLUDE_SUFFIXES):
        return False
    return True


def subtype_order(adata: ad.AnnData) -> list[str]:
    return [str(category) for category in adata.obs[GROUPBY].cat.categories]


def rank_markers(adata: ad.AnnData) -> None:
    sc.tl.rank_genes_groups(
        adata,
        groupby=GROUPBY,
        method="wilcoxon",
        use_raw=True,
        pts=True,
        n_genes=adata.raw.n_vars if adata.raw is not None else adata.n_vars,
        key_added=RANK_KEY,
    )


def group_marker_frame(adata: ad.AnnData, group: str) -> pd.DataFrame:
    rg = adata.uns[RANK_KEY]
    df = pd.DataFrame(
        {
            "gene": pd.Series(rg["names"][group]),
            "score": rg["scores"][group],
            "logfoldchanges": rg["logfoldchanges"][group],
            "pvals": rg["pvals"][group],
            "pvals_adj": rg["pvals_adj"][group],
        }
    )
    pts = rg.get("pts")
    pts_rest = rg.get("pts_rest")
    if pts is not None and pts_rest is not None:
        df["pct_group"] = df["gene"].map(pts[group])
        df["pct_rest"] = df["gene"].map(pts_rest[group])
    else:
        df["pct_group"] = pd.NA
        df["pct_rest"] = pd.NA
    df["pct_diff"] = df["pct_group"] - df["pct_rest"]
    df["specificity_score"] = (
        df["logfoldchanges"].clip(lower=0).fillna(0)
        * df["pct_diff"].clip(lower=0).fillna(0)
        * df["pct_group"].clip(lower=0).fillna(0)
    )
    return df[df["gene"].map(is_informative_gene)].copy()


def select_top_markers(adata: ad.AnnData, groups: list[str]) -> pd.DataFrame:
    selected_rows: list[pd.DataFrame] = []
    used_genes: set[str] = set()

    for group in groups:
        df = group_marker_frame(adata, group)
        strict = df[
            (df["score"] > 0)
            & (df["pvals_adj"] < 0.05)
            & (df["logfoldchanges"] > 0.25)
            & (df["pct_group"] >= 0.10)
            & (df["pct_diff"] >= 0.05)
        ].copy()
        if strict.empty:
            strict = df[
                (df["score"] > 0)
                & (df["logfoldchanges"] > 0.10)
                & (df["pct_group"] >= 0.05)
            ].copy()
        if strict.empty:
            strict = df[df["score"] > 0].copy()

        strict = strict.sort_values(
            ["specificity_score", "logfoldchanges", "score"],
            ascending=[False, False, False],
        )

        group_rows: list[pd.Series] = []
        for _, row in strict.iterrows():
            gene = str(row["gene"])
            if gene in used_genes:
                continue
            group_rows.append(row)
            used_genes.add(gene)
            if len(group_rows) == TOP_N_PER_SUBTYPE:
                break

        if group_rows:
            out = pd.DataFrame(group_rows).copy()
            out.insert(0, "marker_rank", range(1, len(out) + 1))
            out.insert(0, "subtype", group)
            selected_rows.append(out)

    return pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()


def matrix_for_genes(adata: ad.AnnData, genes: list[str]):
    if adata.raw is not None:
        return adata.raw[:, genes].X
    return adata[:, genes].X


def dotplot_stats(adata: ad.AnnData, marker_table: pd.DataFrame, groups: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = marker_table["gene"].astype(str).tolist()
    X = matrix_for_genes(adata, genes)
    group_labels = adata.obs[GROUPBY].astype(str).to_numpy()

    mean_by_group: dict[str, np.ndarray] = {}
    pct_by_group: dict[str, np.ndarray] = {}
    for group in groups:
        mask = group_labels == group
        X_group = X[mask, :]
        if sparse.issparse(X_group):
            mean_by_group[group] = np.asarray(X_group.mean(axis=0)).ravel()
            pct_by_group[group] = np.asarray((X_group > 0).mean(axis=0)).ravel() * 100.0
        else:
            mean_by_group[group] = np.asarray(X_group.mean(axis=0)).ravel()
            pct_by_group[group] = np.asarray((X_group > 0).mean(axis=0)).ravel() * 100.0

    mean_df = pd.DataFrame(mean_by_group, index=genes)
    pct_df = pd.DataFrame(pct_by_group, index=genes)
    scaled = mean_df.copy()
    for gene in genes:
        row = mean_df.loc[gene]
        span = row.max() - row.min()
        scaled.loc[gene] = 0.0 if span == 0 else (row - row.min()) / span
    return scaled, pct_df


def marker_blocks(marker_table: pd.DataFrame) -> list[tuple[str, int, int]]:
    blocks: list[tuple[str, int, int]] = []
    start = 0
    last = str(marker_table.iloc[0]["subtype"])
    for index, subtype in enumerate(marker_table["subtype"].astype(str)):
        if subtype != last:
            blocks.append((last, start, index - 1))
            start = index
            last = subtype
    blocks.append((last, start, len(marker_table) - 1))
    return blocks


def size_from_pct(percent: float) -> float:
    return 8.0 + (percent / 100.0) * 250.0


def save_top_marker_dotplot(adata: ad.AnnData, marker_table: pd.DataFrame, groups: list[str]) -> None:
    scaled, pct_df = dotplot_stats(adata, marker_table, groups)
    genes = marker_table["gene"].astype(str).tolist()
    n_rows = len(genes)
    n_groups = len(groups)

    fig_width = max(13.5, n_groups * 0.70 + 4.8)
    fig_height = max(12.5, n_rows * 0.34 + 3.6)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.14, right=0.74, bottom=0.31, top=0.96)

    x_values: list[int] = []
    y_values: list[int] = []
    colors: list[float] = []
    sizes: list[float] = []
    for row_i, gene in enumerate(genes):
        y = n_rows - 1 - row_i
        for x, group in enumerate(groups):
            x_values.append(x)
            y_values.append(y)
            colors.append(float(scaled.loc[gene, group]))
            sizes.append(size_from_pct(float(pct_df.loc[gene, group])))

    scatter = ax.scatter(
        x_values,
        y_values,
        c=colors,
        s=sizes,
        cmap="Blues",
        vmin=0,
        vmax=1,
        edgecolors="#8c959f",
        linewidths=0.35,
        zorder=2,
    )

    ax.set_xlim(-0.7, n_groups - 0.3)
    ax.set_ylim(-0.7, n_rows - 0.3)
    ax.set_xticks(range(n_groups))
    ax.set_xticklabels(groups, rotation=90, ha="center", va="top")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(list(reversed(genes)))
    ax.set_xlabel("Astrocyte subtype")
    ax.set_ylabel("Top marker gene")
    ax.set_title("Top markers across renamed astrocyte subtypes")
    ax.grid(False)

    legend_percents = [25, 50, 75, 100]
    ax.text(1.06, 0.76, "Percent expressed", transform=ax.transAxes, va="bottom", ha="left", fontsize=10)
    for offset, percent in enumerate(legend_percents):
        y = 0.70 - offset * 0.035
        ax.scatter(
            [1.09],
            [y],
            s=size_from_pct(percent),
            facecolor="#6e7781",
            edgecolor="#6e7781",
            transform=ax.transAxes,
            clip_on=False,
        )
        ax.text(1.125, y, str(percent), transform=ax.transAxes, va="center", ha="left", fontsize=10, rotation=90)

    fig.text(0.80, 0.39, "Average expression", ha="left", va="bottom", fontsize=10)
    cax = fig.add_axes([0.82, 0.18, 0.025, 0.20])
    fig.colorbar(scatter, cax=cax, orientation="vertical")
    fig.savefig(FIG / "astrocyte_refined_c6_c12_merged_top_marker_dotplot.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def get_gene_expr(adata: ad.AnnData, gene: str) -> np.ndarray:
    source = adata.raw.to_adata() if adata.raw is not None else adata
    if gene not in source.var_names:
        raise ValueError(f"{gene} not found in AnnData.var_names")
    x = source[:, [gene]].X
    return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()


def summarize_b2m(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subtype, sub in df.groupby(GROUPBY, observed=True):
        ctrl = sub.loc[sub["group"] == "Ctrl", "B2m"]
        pdv = sub.loc[sub["group"] == "PD", "B2m"]
        if len(ctrl) and len(pdv):
            p_two = stats.mannwhitneyu(pdv, ctrl, alternative="two-sided").pvalue
            p_greater = stats.mannwhitneyu(pdv, ctrl, alternative="greater").pvalue
        else:
            p_two = np.nan
            p_greater = np.nan
        rows.append(
            {
                "astrocyte_subtype": subtype,
                "n_Ctrl": int(len(ctrl)),
                "n_PD": int(len(pdv)),
                "Ctrl_mean_B2m": float(ctrl.mean()) if len(ctrl) else np.nan,
                "PD_mean_B2m": float(pdv.mean()) if len(pdv) else np.nan,
                "PD_minus_Ctrl": float(pdv.mean() - ctrl.mean()) if len(ctrl) and len(pdv) else np.nan,
                "Ctrl_pct_B2m_pos": float((ctrl > 0).mean() * 100) if len(ctrl) else np.nan,
                "PD_pct_B2m_pos": float((pdv > 0).mean() * 100) if len(pdv) else np.nan,
                "cell_level_mwu_p_two_sided": float(p_two) if np.isfinite(p_two) else np.nan,
                "cell_level_mwu_p_PD_greater_exploratory": float(p_greater) if np.isfinite(p_greater) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    valid = out["cell_level_mwu_p_two_sided"].notna()
    out["cell_level_mwu_fdr_two_sided"] = np.nan
    if valid.any():
        out.loc[valid, "cell_level_mwu_fdr_two_sided"] = multipletests(
            out.loc[valid, "cell_level_mwu_p_two_sided"], method="fdr_bh"
        )[1]
    out["stars_two_sided"] = out["cell_level_mwu_p_two_sided"].map(p_to_stars)
    out["stars_PD_greater_exploratory"] = out["cell_level_mwu_p_PD_greater_exploratory"].map(p_to_stars)
    return out


def save_b2m_violin(df: pd.DataFrame, summary: pd.DataFrame, groups: list[str]) -> None:
    fig_w = max(12, 0.82 * len(groups))
    fig, ax = plt.subplots(figsize=(fig_w, 5.8))
    sns.violinplot(
        data=df,
        x=GROUPBY,
        y="B2m",
        hue="group",
        order=groups,
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
        x=GROUPBY,
        y="B2m",
        hue="group",
        order=groups,
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
    for handle, label in zip(handles, labels):
        if label in {"Ctrl", "PD"} and label not in seen:
            keep.append((handle, label))
            seen.add(label)
    ax.legend([h for h, _ in keep], [l for _, l in keep], title="Group", frameon=False, loc="upper right")

    ymax = float(df["B2m"].quantile(0.995))
    y_pad = max(0.12, ymax * 0.08)
    stats_by = summary.set_index("astrocyte_subtype")
    for i, subtype in enumerate(groups):
        if subtype not in stats_by.index:
            continue
        sub = df[df[GROUPBY].astype(str) == subtype]
        local_max = float(sub["B2m"].quantile(0.995)) if len(sub) else ymax
        y = local_max + y_pad
        ax.plot([i - 0.22, i - 0.22, i + 0.22, i + 0.22], [y, y + y_pad * 0.25, y + y_pad * 0.25, y], color="black", lw=0.8)
        ax.text(i, y + y_pad * 0.32, stats_by.loc[subtype, "stars_two_sided"], ha="center", va="bottom", fontsize=9)

    ax.set_ylim(0, max(ymax + y_pad * 3.2, ax.get_ylim()[1]))
    ax.set_xlabel("")
    ax.set_ylabel("B2m expression (log-normalized)")
    ax.set_title("B2m expression across renamed astrocyte subtypes")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        0.01,
        "Stars: cell-level two-sided Mann-Whitney U test; descriptive only because Ctrl/PD each has one library.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIG / "astrocyte_refined_c6_c12_merged_B2m_violin_with_stats.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(H5 / "astrocyte_subclustered_refined.h5ad")
    groups = add_renamed_label(adata)

    rank_markers(adata)
    marker_table = select_top_markers(adata, groups)
    marker_table.to_csv(TAB / "astrocyte_refined_c6_c12_merged_top_markers_for_dotplot.csv", index=False, encoding="utf-8-sig")
    save_top_marker_dotplot(adata, marker_table, groups)

    df = adata.obs[["group", GROUPBY]].copy()
    df["B2m"] = get_gene_expr(adata, "B2m")
    summary = summarize_b2m(df)
    summary["astrocyte_subtype"] = pd.Categorical(summary["astrocyte_subtype"], categories=groups, ordered=True)
    summary = summary.sort_values("astrocyte_subtype")
    summary.to_csv(TAB / "astrocyte_refined_c6_c12_merged_B2m_stats.csv", index=False, encoding="utf-8-sig")
    save_b2m_violin(df, summary, groups)


if __name__ == "__main__":
    main()
