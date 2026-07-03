from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


matplotlib.use("Agg")
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=True)

OUT = Path("p_cresol_snRNA_analysis")
FIG = OUT / "figures"
TAB = OUT / "tables"
H5 = OUT / "h5ad"

GROUPBY = "astrocyte_refined_subtype"
RANK_KEY = "rank_genes_groups_astrocyte_refined_top_markers"
TOP_N_PER_SUBTYPE = 3

EXCLUDE_EXACT = {"Malat1", "Xist"}
EXCLUDE_PREFIXES = ("mt-", "Rpl", "Rps", "Gm")
EXCLUDE_SUFFIXES = ("Rik",)


def is_informative_gene(gene: str) -> bool:
    if gene in EXCLUDE_EXACT:
        return False
    if any(gene.startswith(prefix) for prefix in EXCLUDE_PREFIXES):
        return False
    if any(gene.endswith(suffix) for suffix in EXCLUDE_SUFFIXES):
        return False
    return True


def subtype_order(adata: ad.AnnData) -> list[str]:
    if hasattr(adata.obs[GROUPBY], "cat"):
        return [str(category) for category in adata.obs[GROUPBY].cat.categories]
    return sorted(adata.obs[GROUPBY].astype(str).unique())


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
    names = pd.Series(rg["names"][group], name="gene")
    df = pd.DataFrame(
        {
            "gene": names,
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
    df = df[df["gene"].map(is_informative_gene)].copy()
    return df


def select_top_markers(adata: ad.AnnData, groups: list[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    selected: dict[str, list[str]] = {}
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

        group_genes: list[str] = []
        group_rows: list[pd.Series] = []
        for _, row in strict.iterrows():
            gene = str(row["gene"])
            if gene in used_genes:
                continue
            group_genes.append(gene)
            group_rows.append(row)
            used_genes.add(gene)
            if len(group_genes) == TOP_N_PER_SUBTYPE:
                break

        selected[group] = group_genes
        if group_rows:
            out = pd.DataFrame(group_rows).copy()
            out.insert(0, "marker_rank", range(1, len(out) + 1))
            out.insert(0, "subtype", group)
            selected_rows.append(out)

    table = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()
    return selected, table


def save_dotplot(adata: ad.AnnData, marker_dict: dict[str, list[str]]) -> None:
    marker_dict = {group: genes for group, genes in marker_dict.items() if genes}
    flat_markers = [gene for genes in marker_dict.values() for gene in genes]
    n_genes = sum(len(genes) for genes in marker_dict.values())
    n_groups = len(marker_dict)

    sc.pl.dotplot(
        adata,
        flat_markers,
        groupby=GROUPBY,
        use_raw=True,
        standard_scale="var",
        dendrogram=False,
        swap_axes=True,
        cmap="Blues",
        show=False,
    )
    fig = plt.gcf()
    rotate_scanpy_legend_titles(fig)
    fig.set_size_inches(max(13.5, n_groups * 0.70 + 4.8), max(12.5, n_genes * 0.30 + 4.8))
    fig.savefig(FIG / "astrocyte_refined_top_marker_dotplot.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "astrocyte_refined_top_marker_dotplot_image2_style.png", dpi=300, bbox_inches="tight")
    plt.close("all")


def rotate_scanpy_legend_titles(fig: plt.Figure) -> None:
    for text in fig.findobj(match=matplotlib.text.Text):
        normalized = " ".join(text.get_text().split())
        if normalized == "Fraction of cells in group (%)":
            text.set_rotation(90)
            text.set_ha("center")
            text.set_va("center")
        elif normalized == "Mean expression in group":
            text.set_rotation(90)
            text.set_ha("center")
            text.set_va("center")


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


def save_clear_custom_dotplot(adata: ad.AnnData, marker_table: pd.DataFrame, groups: list[str]) -> None:
    scaled, pct_df = dotplot_stats(adata, marker_table, groups)
    genes = marker_table["gene"].astype(str).tolist()
    n_rows = len(genes)
    n_groups = len(groups)

    fig_width = max(13.5, n_groups * 0.70 + 4.8)
    fig_height = max(12.5, n_rows * 0.34 + 3.6)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.14, right=0.74, bottom=0.34, top=0.98)

    for block_i, (_, start, end) in enumerate(marker_blocks(marker_table)):
        y_top = n_rows - 1 - start + 0.5
        y_bottom = n_rows - 1 - end - 0.5
        if block_i % 2 == 0:
            ax.axhspan(y_bottom, y_top, color="#f7f9fc", zorder=0)
        ax.axhline(y_bottom, color="#d0d7de", linewidth=0.7, zorder=1)

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
    ax.set_xlabel("Astrocyte refined subtype")
    ax.set_ylabel("Top marker gene")
    ax.grid(False)

    legend_percents = [25, 50, 75, 100]
    ax.text(
        1.06,
        0.76,
        "Percent expressed",
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=10,
        clip_on=False,
    )
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
        ax.text(
            1.125,
            y,
            str(percent),
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=10,
            rotation=90,
            clip_on=False,
        )

    fig.text(0.80, 0.39, "Average expression", ha="left", va="bottom", fontsize=10)
    cax = fig.add_axes([0.82, 0.18, 0.025, 0.20])
    cbar = fig.colorbar(scatter, cax=cax, orientation="vertical")

    for filename in [
        "astrocyte_refined_top_marker_dotplot_clear.png",
        "astrocyte_refined_top_marker_dotplot.png",
        "astrocyte_refined_top_marker_dotplot_image2_style.png",
    ]:
        fig.savefig(FIG / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(H5 / "astrocyte_subclustered_refined.h5ad")
    groups = subtype_order(adata)
    rank_markers(adata)
    marker_dict, marker_table = select_top_markers(adata, groups)

    marker_table.to_csv(TAB / "astrocyte_refined_top_markers_for_dotplot.csv", index=False, encoding="utf-8-sig")
    save_dotplot(adata, marker_dict)
    save_clear_custom_dotplot(adata, marker_table, groups)

    (OUT / "ASTRO_TOP_MARKER_DOTPLOT_说明.md").write_text(
        "\n".join(
            [
                "# 星胶 refined 亚群 top marker dotplot",
                "",
                "输出文件：",
                "- `figures/astrocyte_refined_top_marker_dotplot.png`",
                "- `figures/astrocyte_refined_top_marker_dotplot_image2_style.png`",
                "- `figures/astrocyte_refined_top_marker_dotplot_clear.png`",
                "- `tables/astrocyte_refined_top_markers_for_dotplot.csv`",
                "",
                f"方法：对 `astrocyte_refined_subtype` 的每个亚群做 Wilcoxon rank-sum marker 分析，每个亚群优先选取 {TOP_N_PER_SUBTYPE} 个该亚群高表达、其他亚群低表达的 top marker。",
                "",
                "筛选时排除了线粒体、核糖体、Malat1/Xist、Gm* 和 *Rik 等不适合展示为生物学 marker 的基因；dotplot 中点大小表示表达细胞比例，颜色表示亚群间标准化后的平均表达。",
                "",
                "`astrocyte_refined_top_marker_dotplot_clear.png` 为自定义清晰版：横轴使用与 `astrocyte_refined_subtypes_tsne.png` legend 一致的完整亚群名称；右侧不再显示 marker block 文字；图例采用手工排版，`Percent expressed` 横排标题下方竖排圆点和数字，`Average expression` 横排标题下方竖排色柱。",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
