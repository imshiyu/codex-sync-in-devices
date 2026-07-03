from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


matplotlib.use("Agg")
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=True)

OUT = Path("p_cresol_snRNA_analysis")
FIG = OUT / "figures"
TAB = OUT / "tables"
H5 = OUT / "h5ad"

BROAD_MARKERS = {
    "Neuron": ["Snap25", "Syt1", "Rbfox3", "Map2", "Tubb3"],
    "Astrocyte": ["Aldh1l1", "Aqp4", "Slc1a2", "Slc1a3", "Gfap", "S100b"],
    "Microglia": ["C1qa", "C1qb", "C1qc", "P2ry12", "Tmem119", "Aif1"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mobp", "Mag"],
    "OPC": ["Pdgfra", "Cspg4", "Vcan", "Sox10"],
    "Endothelial": ["Pecam1", "Cldn5", "Kdr", "Flt1", "Vwf"],
    "Pericyte_VSMC": ["Pdgfrb", "Rgs5", "Acta2", "Tagln", "Myl9"],
}

NEURON_3_MARKERS = {
    "DA_neuron": ["Th", "Slc6a3", "Ddc", "Slc18a2", "Nr4a2", "Pitx3", "Aldh1a1"],
    "GABAergic": ["Gad1", "Gad2", "Slc32a1", "Pvalb", "Sst"],
    "Glutamatergic": ["Slc17a6", "Slc17a7", "Slc17a8", "Tbr1"],
}


def save_df(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def broad_cell_type(label: str) -> str:
    label = str(label)
    if "neuron" in label.lower():
        return "Neuron"
    if label == "Endothelial_mixed":
        return "Endothelial"
    return label


def set_axes(ax: plt.Axes, basis: str) -> None:
    ax.set_xlabel(f"{basis.upper()}1")
    ax.set_ylabel(f"{basis.upper()}2")
    ax.grid(False)
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=8, length=3)


def plot_embedding(
    adata: ad.AnnData,
    basis: str,
    color: str,
    filename: Path,
    title: str,
    size: float = 4,
    legend_title: str | None = None,
) -> None:
    coords = adata.obsm[f"X_{basis}"]
    obs = adata.obs.copy()
    obs["_x"] = coords[:, 0]
    obs["_y"] = coords[:, 1]
    cats = list(obs[color].astype("category").cat.categories)
    palette = sns.color_palette("tab20", n_colors=max(len(cats), 3))
    lut = dict(zip(cats, palette))

    fig, ax = plt.subplots(figsize=(8.6, 6.6))
    for cat in cats:
        sub = obs[obs[color].astype(str) == str(cat)]
        ax.scatter(sub["_x"], sub["_y"], s=size, c=[lut[cat]], label=str(cat), linewidths=0, alpha=0.82)
    ax.set_title(title)
    set_axes(ax, basis)
    ax.legend(
        title=legend_title or color,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=3,
        fontsize=8,
        title_fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_split_embedding(
    adata: ad.AnnData,
    basis: str,
    color: str,
    split: str,
    filename: Path,
    title: str,
    size: float = 4,
) -> None:
    coords = adata.obsm[f"X_{basis}"]
    obs = adata.obs.copy()
    obs["_x"] = coords[:, 0]
    obs["_y"] = coords[:, 1]
    cats = list(obs[color].astype("category").cat.categories)
    split_vals = list(obs[split].astype("category").cat.categories)
    palette = sns.color_palette("tab20", n_colors=max(len(cats), 3))
    lut = dict(zip(cats, palette))
    fig, axes = plt.subplots(1, len(split_vals), figsize=(6.8 * len(split_vals), 5.8), sharex=True, sharey=True)
    if len(split_vals) == 1:
        axes = [axes]
    for ax, split_val in zip(axes, split_vals):
        sub_obs = obs[obs[split].astype(str) == str(split_val)]
        for cat in cats:
            sub = sub_obs[sub_obs[color].astype(str) == str(cat)]
            ax.scatter(sub["_x"], sub["_y"], s=size, c=[lut[cat]], label=str(cat), linewidths=0, alpha=0.82)
        ax.set_title(str(split_val))
        set_axes(ax, basis)
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, title=color, loc="center left", bbox_to_anchor=(0.995, 0.5), frameon=False, markerscale=3, fontsize=8, title_fontsize=9)
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def cell_type_proportions(adata: ad.AnnData, groupby: str, filename_stem: str) -> pd.DataFrame:
    counts = adata.obs.groupby(["group", groupby], observed=False).size().reset_index(name="n_cells")
    totals = adata.obs.groupby("group", observed=False).size().rename("total_cells").reset_index()
    prop = counts.merge(totals, on="group")
    prop["fraction"] = prop["n_cells"] / prop["total_cells"]
    save_df(prop, TAB / f"{filename_stem}.csv")

    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.barplot(data=prop, x=groupby, y="fraction", hue="group", ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of cells/nuclei")
    ax.tick_params(axis="x", rotation=35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / f"{filename_stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return prop


def dotplot_markers(
    adata: ad.AnnData,
    markers: dict[str, list[str]],
    groupby: str,
    filename: Path,
    title: str,
) -> None:
    present = {
        label: [gene for gene in genes if gene in adata.var_names]
        for label, genes in markers.items()
    }
    present = {label: genes for label, genes in present.items() if genes}
    if not present:
        return
    sc.pl.dotplot(
        adata,
        present,
        groupby=groupby,
        use_raw=True if adata.raw is not None else False,
        standard_scale="var",
        dendrogram=False,
        show=False,
    )
    fig = plt.gcf()
    fig.suptitle(title, y=1.02)
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close("all")


def add_neuron_3subtype(neuron: ad.AnnData) -> pd.Series:
    subtype = neuron.obs["neuron_subtype"].astype(str).copy()
    subtype = subtype.replace({"DA_SNc_Aldh1a1": "DA_neuron", "Neuropeptide": "Other_neuron"})
    subtype = subtype.where(subtype.isin(["DA_neuron", "GABAergic", "Glutamatergic"]), "Other_neuron")
    neuron.obs["neuron_3subtype"] = pd.Categorical(
        subtype,
        categories=["DA_neuron", "GABAergic", "Glutamatergic", "Other_neuron"],
        ordered=False,
    )
    return neuron.obs["neuron_3subtype"].isin(["DA_neuron", "GABAergic", "Glutamatergic"]).to_numpy()


def lightweight_obs_view(adata: ad.AnnData, mask: np.ndarray) -> ad.AnnData:
    obs = adata.obs.loc[mask].copy()
    view = ad.AnnData(X=np.zeros((obs.shape[0], 1), dtype=np.float32), obs=obs)
    for key in ["X_umap", "X_tsne"]:
        if key in adata.obsm:
            view.obsm[key] = np.asarray(adata.obsm[key])[mask]
    return view


def main() -> None:
    main_adata = ad.read_h5ad(H5 / "02_main_annotated.h5ad")
    main_adata.obs["cell_type_original"] = main_adata.obs["cell_type"].astype(str)
    main_adata.obs["cell_type_broad"] = pd.Categorical(main_adata.obs["cell_type_original"].map(broad_cell_type))
    save_df(
        main_adata.obs[["leiden", "cell_type_original", "cell_type_broad"]]
        .drop_duplicates()
        .sort_values(["leiden"]),
        TAB / "major_cluster_annotations_broad.csv",
    )
    cell_type_proportions(main_adata, "cell_type_broad", "cell_proportions_by_cell_type_broad")
    dotplot_markers(
        main_adata,
        BROAD_MARKERS,
        "cell_type_broad",
        FIG / "major_marker_dotplot_broad_cell_type.png",
        "Broad cell type marker dotplot",
    )
    for basis in ["umap", "tsne"]:
        plot_embedding(
            main_adata,
            basis,
            "cell_type_broad",
            FIG / f"main_{basis}_cell_type_broad_axes_legend.png",
            f"Main {basis.upper()} - broad cell types",
            size=4,
            legend_title="Broad cell type",
        )
        plot_embedding(
            main_adata,
            basis,
            "group",
            FIG / f"main_{basis}_group_axes_legend.png",
            f"Main {basis.upper()} - group",
            size=4,
            legend_title="Group",
        )

    main_adata.write_h5ad(H5 / "02_main_annotated_broad.h5ad", compression="gzip")

    neuron = ad.read_h5ad(H5 / "neuron_subclustered.h5ad")
    principal_mask = add_neuron_3subtype(neuron)
    neuron3 = lightweight_obs_view(neuron, principal_mask)
    save_df(
        neuron3.obs[["neuron_leiden", "neuron_subtype", "neuron_3subtype"]]
        .drop_duplicates()
        .sort_values(["neuron_leiden"]),
        TAB / "neuron_subcluster_annotations_3subtype.csv",
    )
    cell_type_proportions(neuron3, "neuron_3subtype", "cell_proportions_by_neuron_3subtype")
    dotplot_markers(
        neuron,
        NEURON_3_MARKERS,
        "neuron_3subtype",
        FIG / "neuron_marker_dotplot_3subtype.png",
        "Neuron three-subtype marker dotplot",
    )
    for basis in ["umap", "tsne"]:
        plot_embedding(
            neuron3,
            basis,
            "neuron_3subtype",
            FIG / f"neuron_{basis}_3subtype_axes_legend.png",
            f"Neuron {basis.upper()} - three principal subtypes",
            size=7,
            legend_title="Neuron subtype",
        )
        plot_split_embedding(
            neuron3,
            basis,
            "neuron_3subtype",
            "group",
            FIG / f"neuron_{basis}_3subtype_by_group_axes_legend.png",
            f"Neuron {basis.upper()} - three principal subtypes by group",
            size=7,
        )
    neuron.write_h5ad(H5 / "neuron_subclustered_3subtype.h5ad", compression="gzip")

    note = OUT / "REVISED_绘图和神经元分群说明.md"
    note.write_text(
        "\n".join(
            [
                "# 绘图和神经元分群修订说明",
                "",
                "按要求已新增 broad cell type 主图：总分群中 neuron 不再细分 DA/GABA/Glut，而统一标为 `Neuron`。",
                "",
                "新增图：",
                "- `figures/main_umap_cell_type_broad_axes_legend.png`",
                "- `figures/main_tsne_cell_type_broad_axes_legend.png`",
                "- `figures/main_umap_group_axes_legend.png`",
                "- `figures/main_tsne_group_axes_legend.png`",
                "- `figures/neuron_umap_3subtype_axes_legend.png`",
                "- `figures/neuron_tsne_3subtype_axes_legend.png`",
                "- `figures/neuron_umap_3subtype_by_group_axes_legend.png`",
                "- `figures/neuron_tsne_3subtype_by_group_axes_legend.png`",
                "- `figures/major_marker_dotplot_broad_cell_type.png`",
                "- `figures/neuron_marker_dotplot_3subtype.png`",
                "",
                "所有新版 UMAP/tSNE 都保留横纵坐标，并把分类文字移到图外 legend，避免遮挡细胞点。",
                "",
                "新增表：",
                "- `tables/major_cluster_annotations_broad.csv`",
                "- `tables/cell_proportions_by_cell_type_broad.csv`",
                "- `tables/neuron_subcluster_annotations_3subtype.csv`",
                "- `tables/cell_proportions_by_neuron_3subtype.csv`",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
