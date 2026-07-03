from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


matplotlib.use("Agg")

OUT = Path("p_cresol_snRNA_analysis")
FIG = OUT / "figures"
TAB = OUT / "tables"
H5 = OUT / "h5ad"

GROUPBY = "astrocyte_refined_subtype"
MERGE_FROM = {
    "AST_B2m_MHCI_IFN_c12",
    "AST_Apoe_Clu_metabolic_c6",
}
MERGED_LABEL = "AST_B2m_MHCI_IFN_1"
RENAMED_COL = "astrocyte_refined_subtype_renamed"

RENAME_MAP = {
    "AST_Aldh1a1_Dao_regional_c1": "AST_regional_1",
    "AST_Aldh1a1_Dao_regional_c2": "AST_regional_2",
    "AST_Aldh1a1_Dao_regional_c3": "AST_regional_3",
    "AST_Aldh1a1_Dao_regional_c5": "AST_regional_4",
    "AST_B2m_MHCI_IFN_c0": "AST_B2m_MHCI_IFN_2",
    "AST_B2m_MHCI_IFN_c10": "AST_B2m_MHCI_IFN_3",
    "AST_Maob_Slc6a11_transport_c11": "AST_transport_1",
    "AST_Maob_Slc6a11_transport_c4": "AST_transport_2",
    "AST_Maob_Slc6a11_transport_c8": "AST_transport_3",
    "AST_glutamate_homeostasis_c7": "AST_glutamate_1",
    "AST_glutamate_homeostasis_c9": "AST_glutamate_2",
}


def set_axes(ax: plt.Axes, basis: str) -> None:
    ax.set_xlabel(f"{basis.upper()}1")
    ax.set_ylabel(f"{basis.upper()}2")
    ax.grid(False)
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=8, length=3)


def renamed_order(original: pd.Series) -> list[str]:
    if hasattr(original, "cat"):
        source_order = [str(category) for category in original.cat.categories]
    else:
        source_order = sorted(original.astype(str).unique())

    order: list[str] = []
    inserted = False
    for label in source_order:
        if label in MERGE_FROM:
            if not inserted:
                order.append(MERGED_LABEL)
                inserted = True
            continue
        order.append(RENAME_MAP.get(label, label))
    return order


def add_renamed_label(adata: ad.AnnData) -> list[str]:
    original = adata.obs[GROUPBY].astype(str)
    renamed = original.map(lambda label: MERGED_LABEL if label in MERGE_FROM else RENAME_MAP.get(label, label))
    order = renamed_order(adata.obs[GROUPBY])
    adata.obs[RENAMED_COL] = pd.Categorical(renamed, categories=order, ordered=True)
    return order


def plot_embedding(adata: ad.AnnData, basis: str, order: list[str]) -> None:
    coords = np.asarray(adata.obsm[f"X_{basis}"])
    obs = adata.obs[[GROUPBY, RENAMED_COL, "group"]].copy()
    obs["_x"] = coords[:, 0]
    obs["_y"] = coords[:, 1]

    palette = sns.color_palette("tab20", n_colors=max(len(order), 3))
    lut = dict(zip(order, palette))

    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    for label in order:
        sub = obs[obs[RENAMED_COL].astype(str) == label]
        if sub.empty:
            continue
        ax.scatter(
            sub["_x"],
            sub["_y"],
            s=8,
            c=[lut[label]],
            label=label,
            linewidths=0,
            alpha=0.78,
            zorder=2,
        )

    ax.set_title(f"Astrocyte refined subtypes: {basis.upper()}")
    set_axes(ax, basis)
    ax.legend(
        title="Astrocyte subtype",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=2.3,
        fontsize=7.3,
        title_fontsize=8.3,
    )
    fig.tight_layout()
    fig.savefig(FIG / f"astrocyte_refined_subtypes_c6_c12_merged_{basis}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_counts(adata: ad.AnnData) -> None:
    counts = (
        adata.obs.groupby([RENAMED_COL, "group"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values([RENAMED_COL, "group"])
    )
    counts.to_csv(TAB / "astrocyte_refined_c6_c12_merged_counts.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(H5 / "astrocyte_subclustered_refined.h5ad")
    order = add_renamed_label(adata)
    save_counts(adata)
    plot_embedding(adata, "umap", order)
    plot_embedding(adata, "tsne", order)


if __name__ == "__main__":
    main()
