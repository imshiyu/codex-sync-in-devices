from __future__ import annotations

import json
import math
import re
import shutil
import zipfile
from collections import OrderedDict
from pathlib import Path

import anndata as ad
import gseapy as gp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


matplotlib.use("Agg")
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)


PROJECT_ROOT = Path("p_cresol_snRNA_analysis")
INPUTS = OrderedDict(
    [
        (
            "Ctrl",
            Path(r"H:\黑质单细胞\data\matrix\25070600_Ctrl\filtered_feature_bc_matrix.zip"),
        ),
        (
            "PD",
            Path(r"H:\黑质单细胞\data\matrix\25070600_PD\filtered_feature_bc_matrix.zip"),
        ),
    ]
)

OUT = PROJECT_ROOT
FIG = OUT / "figures"
TAB = OUT / "tables"
H5 = OUT / "h5ad"
EXTRACT = OUT / "input_10x"
ENRICH = OUT / "enrichment"

RANDOM_STATE = 1337


MARKERS = {
    "Astrocyte": ["Aldh1l1", "Aqp4", "Slc1a2", "Slc1a3", "Gfap", "S100b", "Gja1"],
    "Neuron": ["Snap25", "Syt1", "Rbfox3", "Map2", "Tubb3"],
    "Dopaminergic_neuron": ["Th", "Slc6a3", "Ddc", "Slc18a2", "Nr4a2", "Pitx3", "Aldh1a1"],
    "GABAergic_neuron": ["Gad1", "Gad2", "Slc32a1"],
    "Glutamatergic_neuron": ["Slc17a6", "Slc17a7", "Slc17a8"],
    "Microglia": ["C1qa", "C1qb", "C1qc", "P2ry12", "Tmem119", "Aif1", "Cx3cr1"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mobp", "Mag"],
    "OPC": ["Pdgfra", "Cspg4", "Vcan", "Sox10"],
    "Endothelial": ["Pecam1", "Cldn5", "Kdr", "Flt1", "Vwf"],
    "Pericyte_VSMC": ["Pdgfrb", "Rgs5", "Acta2", "Tagln", "Myl9"],
    "Ependymal": ["Foxj1", "Tmem212", "Dnah5", "Pifo"],
}

ASTRO_PROGRAMS = {
    "B2m_MHCI": ["B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Isg15", "Ifit1"],
    "IFN_ISG": ["Ifnb1", "Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1", "Usp18", "Cxcl10", "Oas1a", "Mx1"],
    "Reactive_complement": ["Gfap", "Lcn2", "Serpina3n", "C3", "C4b", "A2m", "Cd44", "Vim", "Ccl2", "Il6"],
    "Senescence_like": ["Cdkn1a", "Cdkn2a", "Serpine1", "Glb1", "Il6", "Ccl2", "Mmp3", "Lmnb1"],
    "Glutamate_homeostasis": ["Slc1a2", "Slc1a3", "Glul", "Aqp4", "Aldh1l1"],
    "RAGE_ligands": ["B2m", "Hmgb1", "S100b", "S100a6", "S100a8", "S100a9", "S100a10"],
}

NEURON_PROGRAMS = {
    "DA_identity": ["Th", "Slc6a3", "Ddc", "Slc18a2", "Nr4a2", "Pitx3", "Aldh1a1"],
    "Stress_IEG_MAPK": ["Fos", "Jun", "Atf3", "Dusp1", "Dusp5", "Egr1"],
    "NFkB_inflammatory": ["Nfkbia", "Tnfaip3", "Ccl2", "Cxcl10", "Rel", "Rela"],
    "Oxidative_stress": ["Hmox1", "Nqo1", "Sod2", "Gpx4", "Txnrd1", "Gclm"],
    "Apoptosis_DDIT3": ["Bax", "Bbc3", "Pmaip1", "Bcl2l11", "Casp3", "Ddit3"],
}

ASTRO_SUBTYPE_MARKERS = {
    "AST_Apoe_Clu_Cst3_metabolic_homeostatic": ["Apoe", "Clu", "Cst3", "Sparc", "Aldoc"],
    "AST_Aldh1a1_Dao_regional_homeostatic": ["Aldh1a1", "Dao", "Gpc5", "Slc1a3"],
    "AST_Slc1a2_Slc1a3_Glul_glutamate_homeostasis": ["Slc1a2", "Slc1a3", "Glul", "Aqp4", "Aldh1l1"],
    "AST_Gfap_Aqp4_C4b_reactive_complement": ["Gfap", "Aqp4", "C4b", "C3", "Serpina3n", "Lcn2"],
    "AST_B2m_MHCI_IFN_high": ["B2m", "H2-K1", "H2-D1", "Tap1", "Psmb8", "Isg15", "Ifit1"],
    "AST_transport_Maob_Slc6a11": ["Maob", "Slc6a11", "Slc7a10", "Gatm"],
}

NEURON_SUBTYPE_MARKERS = {
    "DA_SNc_Aldh1a1": ["Th", "Slc6a3", "Slc18a2", "Ddc", "Nr4a2", "Pitx3", "Aldh1a1"],
    "GABAergic": ["Gad1", "Gad2", "Slc32a1", "Pvalb", "Sst"],
    "Glutamatergic": ["Slc17a6", "Slc17a7", "Slc17a8", "Tbr1"],
    "Cholinergic": ["Chat", "Slc5a7", "Slc18a3"],
    "Serotonergic": ["Tph2", "Slc6a4", "Ddc"],
    "Neuropeptide": ["Tac1", "Penk", "Pdyn", "Cartpt"],
}

CORE_GENES = OrderedDict(
    [
        ("B2m/MHC-I", ["B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9"]),
        ("IFN/ISG", ["Ifnb1", "Isg15", "Ifit1", "Ifit3", "Cxcl10", "Usp18", "Irf7", "Stat1"]),
        ("Reactive astrocyte", ["Gfap", "Lcn2", "Serpina3n", "C3", "C4b", "A2m", "Cd44", "Vim", "Ccl2", "Il6"]),
        ("Senescence-like", ["Cdkn1a", "Cdkn2a", "Serpine1", "Glb1", "Lmnb1"]),
        ("DA neuron", ["Th", "Slc6a3", "Ddc", "Slc18a2", "Nr4a2", "Pitx3", "Aldh1a1"]),
        ("Neuron stress", ["Nfkbia", "Tnfaip3", "Fos", "Jun", "Atf3", "Hmox1", "Nqo1", "Bax", "Bbc3", "Ddit3"]),
        ("RAGE axis", ["Ager", "B2m", "Hmgb1", "S100b", "S100a6", "S100a8", "S100a9", "S100a10"]),
    ]
)


def ensure_dirs() -> None:
    for d in [OUT, FIG, TAB, H5, EXTRACT, ENRICH]:
        d.mkdir(parents=True, exist_ok=True)


def save_df(df: pd.DataFrame, path: Path, index: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index, encoding="utf-8-sig")


def df_to_md(df: pd.DataFrame, index: bool = False, max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy()
    if index:
        view = view.reset_index()
    columns = [str(c) for c in view.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in view.iterrows():
        vals = []
        for value in row.tolist():
            if isinstance(value, float):
                vals.append(f"{value:.4g}")
            else:
                vals.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_Only first {max_rows} of {len(df)} rows shown._")
    return "\n".join(lines)


def extract_10x_zip(zip_path: Path, sample: str) -> Path:
    dest = EXTRACT / sample
    matrix_file = dest / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    if matrix_file.exists():
        return dest / "filtered_feature_bc_matrix"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest)
    return dest / "filtered_feature_bc_matrix"


def read_inputs() -> ad.AnnData:
    adatas = []
    for group, zip_path in INPUTS.items():
        mtx_dir = extract_10x_zip(zip_path, group)
        a = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True, cache=False)
        a.obs["sample_id"] = f"25070600_{group}"
        a.obs["group"] = group
        a.obs["condition"] = group
        a.obs["barcode"] = a.obs_names.astype(str)
        a.obs_names = [f"{group}_{bc}" for bc in a.obs_names.astype(str)]
        a.var_names_make_unique()
        adatas.append(a)
    combined = ad.concat(adatas, join="outer", label="library", keys=list(INPUTS.keys()), index_unique=None)
    combined.var_names_make_unique()
    return combined


def qc_and_doublets(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.match(r"^Rpl|^Rps", case=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    save_df(adata.obs[["sample_id", "group", "n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]], TAB / "qc_metrics_prefilter.csv")
    qc_summary = adata.obs.groupby("group").agg(
        cells=("group", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
        p95_pct_mt=("pct_counts_mt", lambda x: np.percentile(x, 95)),
    )
    save_df(qc_summary, TAB / "qc_summary_prefilter.csv")

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))
    for ax, y, title in zip(
        axes,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        ["Genes per nucleus", "UMI counts", "Mitochondrial %"],
    ):
        sns.violinplot(data=adata.obs, x="group", y=y, ax=ax, inner="quartile", linewidth=0.7, cut=0)
        ax.set_title(title)
        ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(FIG / "qc_prefilter_violin.png", dpi=300)
    plt.close(fig)

    n_genes_hi = adata.obs.groupby("group")["n_genes_by_counts"].transform(lambda x: np.percentile(x, 99.5))
    counts_hi = adata.obs.groupby("group")["total_counts"].transform(lambda x: np.percentile(x, 99.5))
    qc_mask = (
        (adata.obs["n_genes_by_counts"] >= 200)
        & (adata.obs["total_counts"] >= 500)
        & (adata.obs["pct_counts_mt"] <= 15)
        & (adata.obs["n_genes_by_counts"] <= np.maximum(n_genes_hi, 6000))
        & (adata.obs["total_counts"] <= np.maximum(counts_hi, 30000))
    )
    adata.obs["pass_low_quality_qc"] = qc_mask.values
    adata_qc = adata[qc_mask].copy()

    # Scrublet is run after low-quality filtering, batch-aware by library when available.
    try:
        sc.pp.scrublet(
            adata_qc,
            batch_key="sample_id",
            expected_doublet_rate=0.08,
            sim_doublet_ratio=2.0,
            n_prin_comps=30,
            threshold=0.25,
            random_state=RANDOM_STATE,
        )
    except TypeError:
        sc.pp.scrublet(
            adata_qc,
            expected_doublet_rate=0.08,
            sim_doublet_ratio=2.0,
            n_prin_comps=30,
            threshold=0.25,
            random_state=RANDOM_STATE,
        )
    if "predicted_doublet" not in adata_qc.obs:
        adata_qc.obs["predicted_doublet"] = False
    adata_qc.obs["predicted_doublet"] = adata_qc.obs["predicted_doublet"].fillna(False).astype(bool)

    doublet_summary = adata_qc.obs.groupby("group").agg(
        cells_after_qc=("group", "size"),
        predicted_doublets=("predicted_doublet", "sum"),
        median_doublet_score=("doublet_score", "median"),
    )
    doublet_summary["doublet_rate_pct"] = doublet_summary["predicted_doublets"] / doublet_summary["cells_after_qc"] * 100
    save_df(doublet_summary, TAB / "doublet_summary.csv")

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    sns.histplot(data=adata_qc.obs, x="doublet_score", hue="group", bins=50, element="step", stat="density", common_norm=False, ax=ax)
    ax.set_title("Scrublet doublet scores")
    fig.tight_layout()
    fig.savefig(FIG / "doublet_score_histogram.png", dpi=300)
    plt.close(fig)

    adata_clean = adata_qc[~adata_qc.obs["predicted_doublet"]].copy()
    save_df(
        pd.DataFrame(
            {
                "metric": ["input_cells", "after_low_quality_qc", "after_doublet_filter"],
                **{
                    g: [
                        int((adata.obs["group"] == g).sum()),
                        int((adata_qc.obs["group"] == g).sum()),
                        int((adata_clean.obs["group"] == g).sum()),
                    ]
                    for g in sorted(adata.obs["group"].unique())
                },
            }
        ),
        TAB / "cell_filtering_counts.csv",
        index=False,
    )
    adata_clean.write_h5ad(H5 / "01_qc_doublet_filtered_raw_counts.h5ad", compression="gzip")
    return adata_clean


def normalize_cluster(adata: ad.AnnData, resolution: float = 0.8, prefix: str = "") -> ad.AnnData:
    a = adata.copy()
    if "counts" in a.layers:
        a.X = a.layers["counts"].copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    sc.pp.highly_variable_genes(
        a,
        n_top_genes=min(3000, a.n_vars),
        batch_key="sample_id",
        flavor="seurat",
    )
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, svd_solver="arpack", random_state=RANDOM_STATE, mask_var="highly_variable")
    sc.pp.neighbors(a, n_neighbors=18, n_pcs=35, random_state=RANDOM_STATE)
    sc.tl.umap(a, min_dist=0.35, spread=1.0, random_state=RANDOM_STATE)
    sc.tl.tsne(a, n_pcs=35, perplexity=min(30, max(5, (a.n_obs - 1) // 3)), random_state=RANDOM_STATE)
    sc.tl.leiden(a, resolution=resolution, key_added=f"{prefix}leiden" if prefix else "leiden", random_state=RANDOM_STATE)
    return a


def add_module_scores(adata: ad.AnnData, programs: dict[str, list[str]]) -> None:
    for name, genes in programs.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"{name}_score", random_state=RANDOM_STATE, use_raw=True)
        else:
            adata.obs[f"{name}_score"] = np.nan


def mean_expr_by_group(adata: ad.AnnData, groupby: str, genes: list[str], use_raw: bool = True) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.var_names]
    rows = []
    if not genes:
        return pd.DataFrame()
    source = adata.raw.to_adata() if use_raw and adata.raw is not None else adata
    for group, idx in adata.obs.groupby(groupby).indices.items():
        x = source[idx, genes].X
        if sparse.issparse(x):
            mean_vals = np.asarray(x.mean(axis=0)).ravel()
            pct_vals = np.asarray((x > 0).mean(axis=0)).ravel() * 100
        else:
            mean_vals = x.mean(axis=0)
            pct_vals = (x > 0).mean(axis=0) * 100
        for gene, mean_val, pct_val in zip(genes, mean_vals, pct_vals):
            rows.append({"group": group, "gene": gene, "mean_log_norm": float(mean_val), "pct_expressing": float(pct_val)})
    return pd.DataFrame(rows)


def score_markers_for_clusters(adata: ad.AnnData, cluster_key: str, marker_sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for label, genes in marker_sets.items():
        genes = [g for g in genes if g in adata.var_names]
        if not genes:
            continue
        expr = mean_expr_by_group(adata, cluster_key, genes)
        if expr.empty:
            continue
        scores = expr.groupby("group")["mean_log_norm"].mean()
        pct = expr.groupby("group")["pct_expressing"].mean()
        for cluster in scores.index:
            rows.append(
                {
                    "cluster": str(cluster),
                    "marker_set": label,
                    "score": float(scores.loc[cluster]),
                    "mean_pct_expressing": float(pct.loc[cluster]),
                    "genes_used": ",".join(genes),
                }
            )
    return pd.DataFrame(rows)


def annotate_major_clusters(adata: ad.AnnData, cluster_key: str = "leiden") -> ad.AnnData:
    scores = score_markers_for_clusters(adata, cluster_key, MARKERS)
    save_df(scores, TAB / "major_cluster_marker_scores.csv", index=False)
    pivot = scores.pivot(index="cluster", columns="marker_set", values="score").fillna(0)
    annotations = {}
    for cluster, row in pivot.iterrows():
        ranked = row.sort_values(ascending=False)
        best = ranked.index[0]
        second = ranked.index[1] if len(ranked) > 1 else ""
        if best == "Dopaminergic_neuron":
            label = "Dopaminergic_neuron"
        elif best in {"GABAergic_neuron", "Glutamatergic_neuron"} and row.get("Neuron", 0) > 0.05:
            label = best
        elif best == "Neuron" and row.get("Dopaminergic_neuron", 0) > max(0.05, row.get("GABAergic_neuron", 0), row.get("Glutamatergic_neuron", 0)):
            label = "Dopaminergic_neuron"
        elif best == "Neuron":
            subtype_scores = row[[c for c in ["Dopaminergic_neuron", "GABAergic_neuron", "Glutamatergic_neuron"] if c in row.index]].sort_values(ascending=False)
            label = subtype_scores.index[0] if len(subtype_scores) and subtype_scores.iloc[0] > 0.03 else "Neuron"
        elif best == "OPC" and row.get("Oligodendrocyte", 0) > row.get("OPC", 0) * 0.9:
            label = "Oligodendrocyte_OPC"
        elif best == "Pericyte_VSMC":
            label = "Pericyte_VSMC"
        else:
            label = best
        if ranked.iloc[0] < 0.025:
            label = "Unresolved_low_marker"
        if second and ranked.iloc[1] > ranked.iloc[0] * 0.75 and best not in {"Astrocyte", "Microglia", "Oligodendrocyte", "OPC"}:
            label = f"{label}_mixed"
        annotations[str(cluster)] = label
    adata.obs["cell_type"] = adata.obs[cluster_key].astype(str).map(annotations).astype("category")
    save_df(pd.DataFrame({"cluster": list(annotations), "cell_type": list(annotations.values())}), TAB / "major_cluster_annotations.csv", index=False)
    return adata


def plot_main_embeddings(adata: ad.AnnData, color_keys: list[str], stem: str) -> None:
    for basis in ["umap", "tsne"]:
        sc.pl.embedding(
            adata,
            basis=basis,
            color=color_keys,
            ncols=2,
            size=10 if adata.n_obs < 20000 else 4,
            show=False,
            frameon=False,
        )
        plt.savefig(FIG / f"{stem}_{basis}.png", dpi=300, bbox_inches="tight")
        plt.close("all")


def dotplot_markers(adata: ad.AnnData, marker_dict: dict[str, list[str]], groupby: str, filename: Path) -> None:
    present = OrderedDict((k, [g for g in v if g in adata.var_names]) for k, v in marker_dict.items())
    present = OrderedDict((k, v) for k, v in present.items() if v)
    if not present:
        return
    sc.pl.dotplot(adata, present, groupby=groupby, use_raw=True, standard_scale="var", show=False, dendrogram=False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close("all")


def de_by_group(
    adata: ad.AnnData,
    group_col: str,
    group_a: str = "PD",
    group_b: str = "Ctrl",
    subset_name: str = "all",
    min_pct: float = 0.05,
) -> pd.DataFrame:
    if adata.obs[group_col].nunique() < 2:
        return pd.DataFrame()
    a = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    a.obs = adata.obs.copy()
    # Keep genes expressed in at least a small fraction of the subset.
    x = a.X
    pct_all = np.asarray((x > 0).mean(axis=0)).ravel() if sparse.issparse(x) else (x > 0).mean(axis=0)
    a = a[:, pct_all >= min_pct].copy()
    if a.n_vars < 5:
        return pd.DataFrame()
    sc.tl.rank_genes_groups(a, groupby=group_col, groups=[group_a], reference=group_b, method="wilcoxon", pts=True, use_raw=False)
    df = sc.get.rank_genes_groups_df(a, group=group_a)
    df.insert(0, "subset", subset_name)
    df["direction"] = np.where(df["logfoldchanges"] > 0, "up_in_PD", "down_in_PD")
    return df


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")


def run_enrichr(gene_list: list[str], background_genes: list[str], label: str) -> pd.DataFrame:
    gene_list = sorted({g for g in gene_list if isinstance(g, str) and g})
    if len(gene_list) < 5:
        return pd.DataFrame()
    libraries = ["GO_Biological_Process_2023", "KEGG_2019_Mouse"]
    outdir = ENRICH / safe_name(label)
    outdir.mkdir(parents=True, exist_ok=True)
    results = []
    for lib in libraries:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism="Mouse",
                outdir=str(outdir / safe_name(lib)),
                cutoff=1.0,
                no_plot=True,
                background=background_genes,
            )
            if enr is not None and enr.results is not None and not enr.results.empty:
                df = enr.results.copy()
                df.insert(0, "label", label)
                df.insert(1, "library", lib)
                results.append(df)
        except Exception as exc:
            (outdir / f"{safe_name(lib)}_ERROR.txt").write_text(str(exc), encoding="utf-8")
    if results:
        res = pd.concat(results, ignore_index=True)
        save_df(res, TAB / f"enrichment_{safe_name(label)}.csv", index=False)
        plot_enrichment(res, label)
        return res
    return pd.DataFrame()


def plot_enrichment(enr: pd.DataFrame, label: str, top_n: int = 12) -> None:
    if enr.empty:
        return
    df = enr.copy()
    p_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
    df = df.dropna(subset=[p_col])
    df["neg_log10_p"] = -np.log10(df[p_col].clip(lower=1e-300))
    for lib, sub in df.groupby("library"):
        top = sub.sort_values(p_col).head(top_n).copy()
        if top.empty:
            continue
        top["Term_short"] = top["Term"].str.replace(r"\s*\([^)]*\)$", "", regex=True).str.slice(0, 75)
        top = top.sort_values("neg_log10_p")
        fig, ax = plt.subplots(figsize=(7.2, max(3.2, 0.32 * len(top))))
        ax.barh(top["Term_short"], top["neg_log10_p"], color="#b64b4b" if "up" in label else "#3f6fa8")
        ax.set_xlabel(f"-log10({p_col})")
        ax.set_title(f"{label}: {lib}")
        fig.tight_layout()
        fig.savefig(FIG / f"enrichment_{safe_name(label)}_{safe_name(lib)}.png", dpi=300)
        plt.close(fig)


def cell_type_proportions(adata: ad.AnnData, groupby: str = "cell_type") -> pd.DataFrame:
    counts = adata.obs.groupby(["group", groupby], observed=False).size().reset_index(name="n_cells")
    totals = adata.obs.groupby("group", observed=False).size().rename("total_cells").reset_index()
    prop = counts.merge(totals, on="group")
    prop["fraction"] = prop["n_cells"] / prop["total_cells"]
    save_df(prop, TAB / f"cell_proportions_by_{groupby}.csv", index=False)
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    plot_df = prop.copy()
    plot_df[groupby] = plot_df[groupby].astype(str)
    sns.barplot(data=plot_df, x=groupby, y="fraction", hue="group", ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("Fraction of cells/nuclei")
    ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(FIG / f"cell_proportions_by_{groupby}.png", dpi=300)
    plt.close(fig)
    return prop


def compare_gene_or_score_by_group(adata: ad.AnnData, features: list[str], label: str, obs_group: str = "group") -> pd.DataFrame:
    rows = []
    source = adata.raw.to_adata() if adata.raw is not None else adata
    for feat in features:
        if feat in adata.obs.columns:
            vals = adata.obs[feat].astype(float)
            for group, sub in adata.obs.groupby(obs_group, observed=False):
                rows.append({"feature": feat, "group": group, "mean": float(vals.loc[sub.index].mean()), "median": float(vals.loc[sub.index].median()), "pct_positive": np.nan})
            if set(adata.obs[obs_group].unique()) >= {"Ctrl", "PD"}:
                ctrl = vals[adata.obs[obs_group] == "Ctrl"]
                pdv = vals[adata.obs[obs_group] == "PD"]
                _, p = stats.mannwhitneyu(pdv, ctrl, alternative="two-sided")
                rows.append({"feature": feat, "group": "PD_vs_Ctrl_p", "mean": float(p), "median": np.nan, "pct_positive": np.nan})
        elif feat in adata.var_names:
            x = source[:, [feat]].X
            arr = np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()
            vals = pd.Series(arr, index=adata.obs_names)
            for group, sub in adata.obs.groupby(obs_group, observed=False):
                v = vals.loc[sub.index]
                rows.append({"feature": feat, "group": group, "mean": float(v.mean()), "median": float(v.median()), "pct_positive": float((v > 0).mean() * 100)})
            if set(adata.obs[obs_group].unique()) >= {"Ctrl", "PD"}:
                ctrl = vals[adata.obs[obs_group] == "Ctrl"]
                pdv = vals[adata.obs[obs_group] == "PD"]
                _, p = stats.mannwhitneyu(pdv, ctrl, alternative="two-sided")
                rows.append({"feature": feat, "group": "PD_vs_Ctrl_p", "mean": float(p), "median": np.nan, "pct_positive": np.nan})
    df = pd.DataFrame(rows)
    save_df(df, TAB / f"{safe_name(label)}_group_comparison.csv", index=False)
    return df


def subset_recluster(
    adata: ad.AnnData,
    mask: pd.Series | np.ndarray,
    name: str,
    resolution: float,
    marker_sets: dict[str, list[str]],
    programs: dict[str, list[str]] | None = None,
) -> ad.AnnData | None:
    if int(np.sum(mask)) < 80:
        return None
    sub = adata[mask].copy()
    sub = normalize_cluster(sub, resolution=resolution, prefix=f"{safe_name(name)}_")
    cluster_key = f"{safe_name(name)}_leiden"
    add_module_scores(sub, programs or {})
    scores = score_markers_for_clusters(sub, cluster_key, marker_sets)
    save_df(scores, TAB / f"{safe_name(name)}_subcluster_marker_scores.csv", index=False)
    pivot = scores.pivot(index="cluster", columns="marker_set", values="score").fillna(0)
    annotations = {}
    for cluster, row in pivot.iterrows():
        best = row.sort_values(ascending=False).index[0] if len(row) else "Unresolved"
        annotations[str(cluster)] = best if row.max() >= 0.02 else "Unresolved"
    sub.obs[f"{safe_name(name)}_subtype"] = sub.obs[cluster_key].astype(str).map(annotations).astype("category")
    save_df(pd.DataFrame({"cluster": list(annotations), "subtype": list(annotations.values())}), TAB / f"{safe_name(name)}_subcluster_annotations.csv", index=False)
    plot_main_embeddings(sub, [cluster_key, f"{safe_name(name)}_subtype", "group"], f"{safe_name(name)}_subclusters")
    dotplot_markers(sub, marker_sets, f"{safe_name(name)}_subtype", FIG / f"{safe_name(name)}_subcluster_marker_dotplot.png")
    cell_type_proportions(sub, groupby=f"{safe_name(name)}_subtype")
    sub.write_h5ad(H5 / f"{safe_name(name)}_subclustered.h5ad", compression="gzip")
    return sub


def astrocyte_deep_dive(adata: ad.AnnData) -> ad.AnnData | None:
    astro_mask = adata.obs["cell_type"].astype(str).str.contains("Astrocyte", case=False, na=False)
    astro = subset_recluster(adata, astro_mask, "astrocyte", 1.2, ASTRO_SUBTYPE_MARKERS, ASTRO_PROGRAMS)
    if astro is None:
        return None
    subtype_key = "astrocyte_subtype"

    # Summarize B2m and project-specific modules globally and by astrocyte subtype.
    features = ["B2m", "Ager"] + [f"{k}_score" for k in ASTRO_PROGRAMS]
    compare_gene_or_score_by_group(astro, features, "astrocyte_overall")
    rows = []
    for subtype in astro.obs[subtype_key].cat.categories:
        sub = astro[astro.obs[subtype_key] == subtype].copy()
        comp = compare_gene_or_score_by_group(sub, features, f"astrocyte_{subtype}")
        comp.insert(0, "astro_subtype", subtype)
        rows.append(comp)
    if rows:
        all_comp = pd.concat(rows, ignore_index=True)
        save_df(all_comp, TAB / "astrocyte_subtype_gene_score_group_comparison.csv", index=False)

    b2m_table = []
    source = astro.raw.to_adata() if astro.raw is not None else astro
    if "B2m" in astro.var_names:
        vals = np.asarray(source[:, ["B2m"]].X.toarray()).ravel() if sparse.issparse(source[:, ["B2m"]].X) else np.asarray(source[:, ["B2m"]].X).ravel()
        astro.obs["B2m_expr"] = vals
        for subtype, df in astro.obs.groupby(subtype_key, observed=False):
            ctrl = df.loc[df["group"] == "Ctrl", "B2m_expr"]
            pdv = df.loc[df["group"] == "PD", "B2m_expr"]
            if len(ctrl) and len(pdv):
                stat, p = stats.mannwhitneyu(pdv, ctrl, alternative="two-sided")
                b2m_table.append(
                    {
                        "astro_subtype": subtype,
                        "Ctrl_mean_B2m": float(ctrl.mean()),
                        "PD_mean_B2m": float(pdv.mean()),
                        "PD_minus_Ctrl": float(pdv.mean() - ctrl.mean()),
                        "Ctrl_pct_B2m_pos": float((ctrl > 0).mean() * 100),
                        "PD_pct_B2m_pos": float((pdv > 0).mean() * 100),
                        "n_Ctrl": int(len(ctrl)),
                        "n_PD": int(len(pdv)),
                        "cell_level_mwu_p_descriptive": float(p),
                    }
                )
    b2m_df = pd.DataFrame(b2m_table)
    if not b2m_df.empty:
        b2m_df["cell_level_mwu_fdr_descriptive"] = multipletests(b2m_df["cell_level_mwu_p_descriptive"], method="fdr_bh")[1]
        b2m_df = b2m_df.sort_values(["PD_minus_Ctrl", "PD_mean_B2m"], ascending=False)
        save_df(b2m_df, TAB / "astrocyte_subtype_B2m_PD_vs_Ctrl.csv", index=False)
        fig, ax = plt.subplots(figsize=(9, max(3.5, 0.35 * len(b2m_df))))
        sns.barplot(data=b2m_df, y="astro_subtype", x="PD_minus_Ctrl", ax=ax, color="#9b3f3f")
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("PD - Ctrl mean B2m (log-normalized)")
        ax.set_ylabel("")
        fig.tight_layout()
        fig.savefig(FIG / "astrocyte_subtype_B2m_PD_minus_Ctrl.png", dpi=300)
        plt.close(fig)

    astro_de = de_by_group(astro, "group", subset_name="astrocyte_all")
    save_df(astro_de, TAB / "DE_astrocyte_PD_vs_Ctrl.csv", index=False)
    background = list(astro.var_names)
    up = astro_de.query("logfoldchanges > 0 and pvals_adj < 0.1")["names"].tolist() if not astro_de.empty else []
    down = astro_de.query("logfoldchanges < 0 and pvals_adj < 0.1")["names"].tolist() if not astro_de.empty else []
    if len(up) < 10 and not astro_de.empty:
        up = astro_de.sort_values("scores", ascending=False).head(200)["names"].tolist()
    if len(down) < 10 and not astro_de.empty:
        down = astro_de.sort_values("scores", ascending=True).head(200)["names"].tolist()
    run_enrichr(up, background, "astrocyte_up_in_PD")
    run_enrichr(down, background, "astrocyte_down_in_PD")
    return astro


def neuron_deep_dive(adata: ad.AnnData) -> ad.AnnData | None:
    neuron_mask = adata.obs["cell_type"].astype(str).str.contains("neuron", case=False, na=False)
    neuron = subset_recluster(adata, neuron_mask, "neuron", 1.0, NEURON_SUBTYPE_MARKERS, NEURON_PROGRAMS)
    if neuron is None:
        return None
    subtype_key = "neuron_subtype"
    features = [f"{k}_score" for k in NEURON_PROGRAMS] + ["Th", "Slc6a3", "Slc18a2", "Nr4a2", "Pitx3", "Aldh1a1", "B2m", "Ager"]
    compare_gene_or_score_by_group(neuron, features, "neuron_overall")
    rows = []
    for subtype in neuron.obs[subtype_key].cat.categories:
        sub = neuron[neuron.obs[subtype_key] == subtype].copy()
        comp = compare_gene_or_score_by_group(sub, features, f"neuron_{subtype}")
        comp.insert(0, "neuron_subtype", subtype)
        rows.append(comp)
        de = de_by_group(sub, "group", subset_name=f"neuron_{subtype}")
        save_df(de, TAB / f"DE_neuron_{safe_name(subtype)}_PD_vs_Ctrl.csv", index=False)
        background = list(sub.var_names)
        if not de.empty:
            up = de.query("logfoldchanges > 0 and pvals_adj < 0.1")["names"].tolist()
            down = de.query("logfoldchanges < 0 and pvals_adj < 0.1")["names"].tolist()
            if len(up) < 10:
                up = de.sort_values("scores", ascending=False).head(150)["names"].tolist()
            if len(down) < 10:
                down = de.sort_values("scores", ascending=True).head(150)["names"].tolist()
            run_enrichr(up, background, f"neuron_{subtype}_up_in_PD")
            run_enrichr(down, background, f"neuron_{subtype}_down_in_PD")
    if rows:
        save_df(pd.concat(rows, ignore_index=True), TAB / "neuron_subtype_gene_score_group_comparison.csv", index=False)
    de_all = de_by_group(neuron, "group", subset_name="neuron_all")
    save_df(de_all, TAB / "DE_neuron_all_PD_vs_Ctrl.csv", index=False)
    return neuron


def major_de_and_enrichment(adata: ad.AnnData) -> None:
    background = list(adata.var_names)
    all_enr = []
    for cell_type in sorted(adata.obs["cell_type"].astype(str).unique()):
        mask = adata.obs["cell_type"].astype(str) == cell_type
        if int(mask.sum()) < 50:
            continue
        sub = adata[mask].copy()
        de = de_by_group(sub, "group", subset_name=cell_type)
        save_df(de, TAB / f"DE_{safe_name(cell_type)}_PD_vs_Ctrl.csv", index=False)
        if de.empty:
            continue
        up = de.query("logfoldchanges > 0 and pvals_adj < 0.1")["names"].tolist()
        down = de.query("logfoldchanges < 0 and pvals_adj < 0.1")["names"].tolist()
        if len(up) < 10:
            up = de.sort_values("scores", ascending=False).head(200)["names"].tolist()
        if len(down) < 10:
            down = de.sort_values("scores", ascending=True).head(200)["names"].tolist()
        for direction, genes in [("up_in_PD", up), ("down_in_PD", down)]:
            enr = run_enrichr(genes, background, f"{cell_type}_{direction}")
            if not enr.empty:
                all_enr.append(enr)
    if all_enr:
        save_df(pd.concat(all_enr, ignore_index=True), TAB / "enrichment_all_major_celltypes.csv", index=False)


def core_gene_heatmap(adata: ad.AnnData) -> None:
    genes = []
    for group_genes in CORE_GENES.values():
        genes.extend(group_genes)
    genes = [g for g in OrderedDict.fromkeys(genes) if g in adata.var_names]
    if not genes:
        return
    mean_df = mean_expr_by_group(adata, "cell_type", genes, use_raw=True)
    save_df(mean_df, TAB / "core_gene_mean_expression_by_celltype.csv", index=False)
    mat = mean_df.pivot(index="gene", columns="group", values="mean_log_norm").reindex(genes)
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(max(6, 0.34 * z.shape[1]), max(8, 0.22 * z.shape[0])))
    sns.heatmap(z, cmap="vlag", center=0, ax=ax, cbar_kws={"label": "row z-score"})
    ax.set_title("Core project genes by annotated cell type")
    fig.tight_layout()
    fig.savefig(FIG / "core_gene_heatmap_by_celltype.png", dpi=300)
    plt.close(fig)


def write_summary(adata: ad.AnnData, astro: ad.AnnData | None, neuron: ad.AnnData | None) -> None:
    filter_counts = pd.read_csv(TAB / "cell_filtering_counts.csv")
    props = pd.read_csv(TAB / "cell_proportions_by_cell_type.csv") if (TAB / "cell_proportions_by_cell_type.csv").exists() else pd.DataFrame()
    astro_b2m_path = TAB / "astrocyte_subtype_B2m_PD_vs_Ctrl.csv"
    astro_b2m = pd.read_csv(astro_b2m_path) if astro_b2m_path.exists() else pd.DataFrame()
    major_ann = pd.read_csv(TAB / "major_cluster_annotations.csv") if (TAB / "major_cluster_annotations.csv").exists() else pd.DataFrame()
    qc_summary = pd.read_csv(TAB / "qc_summary_prefilter.csv")
    doublets = pd.read_csv(TAB / "doublet_summary.csv")

    lines = []
    lines.append("# p-cresol PD 黑质 snRNA-seq 分析报告\n")
    lines.append("## 数据和流程\n")
    lines.append("- 输入：`H:\\黑质单细胞\\data\\matrix\\25070600_Ctrl\\filtered_feature_bc_matrix.zip` 与 `25070600_PD` filtered 10X matrix。")
    lines.append("- 工具：Scanpy/AnnData；Scrublet 用于双细胞预测；Leiden 聚类；UMAP 和 tSNE 降维。")
    lines.append("- 过滤：`n_genes_by_counts >= 200`、`total_counts >= 500`、`pct_counts_mt <= 15`，并去除每组极端高基因/UMI细胞和 Scrublet 双细胞。")
    lines.append("- 统计限制：当前目录中 Ctrl 与 PD 各一个 10X library，没有生物学重复；细胞级 p 值仅作描述性排序，不能替代动物/样本层面的显著性检验。\n")

    lines.append("## QC 概览\n")
    lines.append(df_to_md(filter_counts, index=False))
    lines.append("\n预过滤 QC：\n")
    lines.append(df_to_md(qc_summary, index=False))
    lines.append("\n双细胞预测：\n")
    lines.append(df_to_md(doublets, index=False))

    lines.append("\n## 主要细胞类型注释\n")
    if not major_ann.empty:
        lines.append(df_to_md(major_ann, index=False))
    lines.append("\n主要图：")
    lines.append("- `figures/main_umap.png`")
    lines.append("- `figures/main_tsne.png`")
    lines.append("- `figures/major_marker_dotplot.png`")
    lines.append("- `figures/cell_proportions_by_cell_type.png`\n")

    if not props.empty:
        lines.append("细胞比例表见 `tables/cell_proportions_by_cell_type.csv`。由于每组只有一个样本，此处解释为组成变化描述。")

    lines.append("\n## 星形胶质细胞重点结果\n")
    if astro is not None:
        lines.append("- 输出：`h5ad/astrocyte_subclustered.h5ad`、`figures/astrocyte_subclusters_umap.png`、`figures/astrocyte_subcluster_marker_dotplot.png`。")
        if not astro_b2m.empty:
            lines.append("- 星胶亚群 B2m 的 PD-Ctrl 差值如下，优先关注真实出现 PD 上升的亚群：\n")
            cols = ["astro_subtype", "Ctrl_mean_B2m", "PD_mean_B2m", "PD_minus_Ctrl", "Ctrl_pct_B2m_pos", "PD_pct_B2m_pos", "n_Ctrl", "n_PD"]
            lines.append(df_to_md(astro_b2m[cols], index=False))
            up = astro_b2m[astro_b2m["PD_minus_Ctrl"] > 0]
            if not up.empty:
                best = up.iloc[0]
                lines.append(f"\n- B2m 在 `{best['astro_subtype']}` 亚群中 PD 高于 Ctrl；这比强行给全体星胶下结论更稳妥。")
            else:
                lines.append("\n- 本轮星胶亚群中未看到 B2m 的 PD-Ctrl 平均表达差值为正；不建议写 PD 星胶 B2m 上升。")
        lines.append("- 星胶项目基因/score 统计：`tables/astrocyte_subtype_gene_score_group_comparison.csv`。")
        lines.append("- 星胶 GO/KEGG：`tables/enrichment_astrocyte_up_in_PD.csv` 与 `tables/enrichment_astrocyte_down_in_PD.csv`。")
    else:
        lines.append("- 未获得足够星胶细胞用于重聚类。")

    lines.append("\n## 神经元重点结果\n")
    if neuron is not None:
        lines.append("- 输出：`h5ad/neuron_subclustered.h5ad`、`figures/neuron_subclusters_umap.png`、`figures/neuron_subcluster_marker_dotplot.png`。")
        lines.append("- DA identity、stress、NF-kB、oxidative stress、apoptosis 分数见 `tables/neuron_subtype_gene_score_group_comparison.csv`。")
        lines.append("- 各神经元亚群差异和 GO/KEGG 富集见 `tables/DE_neuron_*` 与 `tables/enrichment_neuron_*`。")
    else:
        lines.append("- 未获得足够神经元用于重聚类。")

    lines.append("\n## RAGE/Ager 解释边界\n")
    lines.append("本流程会输出 `Ager` 在各细胞和星胶/神经元亚群中的表达。若 `Ager` 低或未检出，应按项目历史中的保守表述：转录组不能单独证明 RAGE 受体端变化，需要 RAGE IF/WB/RNAscope 或阻断实验。")

    lines.append("\n## 文件索引\n")
    lines.append("- 过滤后主对象：`h5ad/02_main_annotated.h5ad`")
    lines.append("- QC 表：`tables/qc_summary_prefilter.csv`、`tables/doublet_summary.csv`")
    lines.append("- 注释表：`tables/major_cluster_annotations.csv`、`tables/major_cluster_marker_scores.csv`")
    lines.append("- 差异表：`tables/DE_*_PD_vs_Ctrl.csv`")
    lines.append("- 富集表：`tables/enrichment_*.csv`")
    lines.append("- 图片：`figures/*.png`")

    (OUT / "p_cresol_snRNA_analysis_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    adata = read_inputs()
    adata_clean = qc_and_doublets(adata)
    adata_main = normalize_cluster(adata_clean, resolution=0.8)
    add_module_scores(adata_main, {**ASTRO_PROGRAMS, **NEURON_PROGRAMS})
    annotate_major_clusters(adata_main)
    plot_main_embeddings(adata_main, ["leiden", "cell_type", "group"], "main")
    dotplot_markers(adata_main, MARKERS, "leiden", FIG / "major_marker_dotplot_by_cluster.png")
    dotplot_markers(adata_main, MARKERS, "cell_type", FIG / "major_marker_dotplot.png")
    cell_type_proportions(adata_main, "cell_type")
    core_gene_heatmap(adata_main)
    major_de_and_enrichment(adata_main)
    astro = astrocyte_deep_dive(adata_main)
    neuron = neuron_deep_dive(adata_main)
    adata_main.write_h5ad(H5 / "02_main_annotated.h5ad", compression="gzip")
    meta = {
        "random_state": RANDOM_STATE,
        "scanpy_version": sc.__version__,
        "inputs": {k: str(v) for k, v in INPUTS.items()},
        "notes": [
            "Ctrl and PD each have one library; cell-level tests are descriptive.",
            "Astrocyte B2m conclusions should be drawn from observed global/subtype summaries, not forced annotation.",
        ],
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    write_summary(adata_main, astro, neuron)


if __name__ == "__main__":
    main()
