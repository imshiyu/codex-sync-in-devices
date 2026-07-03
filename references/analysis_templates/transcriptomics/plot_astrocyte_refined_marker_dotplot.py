from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc


matplotlib.use("Agg")
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=True)

OUT = Path("p_cresol_snRNA_analysis")
FIG = OUT / "figures"
H5 = OUT / "h5ad"

ASTRO_MARKERS = {
    "Pan-astrocyte": ["Aldh1l1", "Aqp4", "Slc1a2", "Slc1a3", "Gja1", "S100b"],
    "B2m/MHC-I/IFN": ["B2m", "H2-K1", "H2-D1", "Tap1", "Psmb8", "Psmb9", "Isg15", "Ifit1", "Stat1"],
    "Reactive/complement": ["Gfap", "C3", "C4b", "Serpina3n", "Lcn2", "A2m", "Cd44", "Vim"],
    "Glutamate homeostasis": ["Slc1a2", "Slc1a3", "Glul", "Aqp4"],
    "Regional Aldh1a1/Dao": ["Aldh1a1", "Dao", "Gpc5"],
    "Transport Maob/Slc6a11": ["Maob", "Slc6a11", "Slc7a10", "Gatm"],
    "Apoe/Clu metabolic": ["Apoe", "Clu", "Cst3", "Sparc", "Aldoc"],
    "Contamination check": ["Mbp", "Plp1", "Pdgfra", "C1qa", "Pecam1", "Snap25"],
}

SUBTYPE_SHORT_LABELS = {
    "AST_Aldh1a1_Dao_regional_c1": "c1 Aldh1a1/Dao",
    "AST_Aldh1a1_Dao_regional_c2": "c2 Aldh1a1/Dao",
    "AST_Aldh1a1_Dao_regional_c3": "c3 Aldh1a1/Dao",
    "AST_Aldh1a1_Dao_regional_c5": "c5 Aldh1a1/Dao",
    "AST_Apoe_Clu_metabolic_c6": "c6 Apoe/Clu",
    "AST_B2m_MHCI_IFN_c0": "c0 B2m/MHC-I",
    "AST_B2m_MHCI_IFN_c10": "c10 B2m/MHC-I",
    "AST_B2m_MHCI_IFN_c12": "c12 B2m/MHC-I",
    "AST_Maob_Slc6a11_transport_c11": "c11 Maob/Slc6a11",
    "AST_Maob_Slc6a11_transport_c4": "c4 Maob/Slc6a11",
    "AST_Maob_Slc6a11_transport_c8": "c8 Maob/Slc6a11",
    "AST_glutamate_homeostasis_c7": "c7 Glutamate",
    "AST_glutamate_homeostasis_c9": "c9 Glutamate",
}


def flatten_markers(adata: ad.AnnData) -> list[str]:
    gene_index = adata.raw.var_names if adata.raw is not None else adata.var_names
    flat_markers: list[str] = []
    seen: set[str] = set()
    for genes in ASTRO_MARKERS.values():
        for gene in genes:
            if gene in gene_index and gene not in seen:
                flat_markers.append(gene)
                seen.add(gene)
    return flat_markers


def add_short_subtype_labels(adata: ad.AnnData) -> str:
    source = "astrocyte_refined_subtype"
    target = "astrocyte_refined_subtype_short"
    labels = adata.obs[source].map(SUBTYPE_SHORT_LABELS).fillna(adata.obs[source].astype(str))
    if hasattr(adata.obs[source], "cat"):
        ordered = [
            SUBTYPE_SHORT_LABELS.get(str(category), str(category))
            for category in adata.obs[source].cat.categories
        ]
    else:
        ordered = sorted(labels.unique())
    adata.obs[target] = labels
    adata.obs[target] = adata.obs[target].astype("category")
    adata.obs[target] = adata.obs[target].cat.reorder_categories(ordered, ordered=True)
    return target


def main() -> None:
    adata = ad.read_h5ad(H5 / "astrocyte_subclustered_refined.h5ad")
    markers = flatten_markers(adata)
    short_groupby = add_short_subtype_labels(adata)
    sc.pl.dotplot(
        adata,
        markers,
        groupby="astrocyte_refined_subtype",
        use_raw=True,
        standard_scale="var",
        dendrogram=False,
        swap_axes=False,
        show=False,
    )
    fig = plt.gcf()
    fig.set_size_inches(15.5, 6.8)
    fig.savefig(FIG / "astrocyte_refined_cell_type_marker_dotplot.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "all-dotplot.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    sc.pl.dotplot(
        adata,
        markers,
        groupby=short_groupby,
        use_raw=True,
        standard_scale="var",
        dendrogram=False,
        swap_axes=True,
        show=False,
    )
    fig = plt.gcf()
    fig.set_size_inches(9.5, 13.5)
    fig.savefig(FIG / "astrocyte_refined_marker_dotplot_image2_style.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    (OUT / "ASTRO_MARKER_DOTPLOT_说明.md").write_text(
        "\n".join(
            [
                "# 星胶 refined 亚群 marker dotplot",
                "",
                "新增输出：",
                "- `figures/astrocyte_refined_cell_type_marker_dotplot.png`",
                "- `figures/all-dotplot.png`",
                "- `figures/astrocyte_refined_marker_dotplot_image2_style.png`",
                "",
                "本版按扁平基因列表绘制，去掉 dotplot 顶部的大类分组标题。基因顺序保留原 marker 体系顺序，并对重复基因去重。",
                "",
                "`astrocyte_refined_marker_dotplot_image2_style.png` 为参考 `图片2.png` 的竖版 dotplot：Features 在纵轴，astrocyte refined subtype 在横轴，并使用短标签显示具体分群。",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
