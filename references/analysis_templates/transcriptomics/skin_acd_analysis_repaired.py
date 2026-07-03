#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function

"""
Repaired analysis entrypoint for the OXAURU_FC1.5 project.

This script preserves the copied results in `Results/` and writes any newly
generated outputs into `Results_rebuilt/` by default.
"""

import argparse
import os
import shutil
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle
from scipy.stats import ttest_ind


plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


HEATMAP_COLORS = [
    "#D73027", "#F46D43", "#FDAE61", "#FEE090",
    "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4",
]
HEATMAP_CMAP = LinearSegmentedColormap.from_list("custom_rdbu", HEATMAP_COLORS)
HEATMAP_CMAP_R = HEATMAP_CMAP.reversed()

DATASET_COLORS = {"DMF": "#E74C3C", "OXA": "#3498DB", "URU": "#2ECC71"}

TH_COLORS = {
    "Th1": "#E74C3C",
    "Th2": "#3498DB",
    "Th17": "#F39C12",
    "Th22": "#9B59B6",
    "Treg": "#2ECC71",
    "Tfh": "#E67E22",
    "Th9": "#1ABC9C",
    "Tr1": "#95A5A6",
}

TH_MARKERS_MOUSE = {
    "Th1": [
        "Tbx21", "Ifng", "Cxcr3", "Ccr5", "Il12rb2", "Stat4", "Il18r1",
        "Faslg", "Tnf", "Lta", "Il2", "Stat1", "Irf1", "Cxcl9", "Cxcl10",
        "Cxcl11", "Gzmb", "Prf1", "Mx1",
    ],
    "Th2": [
        "Gata3", "Il4", "Il5", "Il13", "Ptgdr2", "Il1rl1", "Ccr4", "Ccr3",
        "Ccr8", "Il17rb", "Stat6", "Areg", "Il25", "Il31", "Tslp", "Ccl17",
        "Ccl22", "Hpgds", "Ccl5",
    ],
    "Th17": ["Il17f", "Il23a", "Ccl20", "Lcn2", "Cxcl1", "Cxcl2", "S100a9"],
    "Th22": ["Ahr", "Il22", "Ccr10", "Ccr6", "Ccr4", "Fgf2", "Il13", "Tnf"],
    "Treg": [
        "Foxp3", "Il2ra", "Ctla4", "Entpd1", "Nt5e", "Ikzf2", "Tnfrsf18",
        "Nrp1", "Tigit", "Lag3", "Il10", "Tgfb1", "Ebi3", "Lrrc32",
    ],
    "Tfh": [
        "Bcl6", "Cxcr5", "Cxcl13", "Il21", "Ascl2", "Pdcd1", "Icos",
        "Sh2d1a", "Btla", "Maf", "Tox2", "Il4", "Cxcr4",
    ],
    "Th9": ["Spi1", "Il9", "Irf4", "Batf", "Il10", "Foxo1", "Smad3"],
    "Tr1": ["Maf", "Il10", "Lag3", "Itga2", "Prdm1", "Ahr", "Havcr2", "Eomes"],
}

CTRL_SAMPLES = ["GSM3832853", "GSM3832854", "GSM3832855"]
OXA_SAMPLES = ["GSM3832856", "GSM3832857", "GSM3832858"]
URU_SAMPLES = ["GSM3832859", "GSM3832860", "GSM3832861"]


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def require_columns(df, required, label):
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError("{0} is missing required columns: {1}".format(
            label, ", ".join(missing)
        ))


def clean_symbol(value):
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def normalize_symbol_series(series):
    return series.fillna("").astype(str).map(clean_symbol)


def deduplicate_by_abs_fc(df, symbol_col, fc_col):
    working = df.copy()
    working[symbol_col] = normalize_symbol_series(working[symbol_col])
    working = working[working[symbol_col] != ""].copy()
    working["_abs_fc"] = working[fc_col].abs()
    working = working.sort_values("_abs_fc", ascending=False)
    working = working.drop_duplicates(symbol_col, keep="first")
    return working.drop("_abs_fc", axis=1)


def build_symbol_to_fc(df, symbol_col, fc_col):
    unique_df = deduplicate_by_abs_fc(df, symbol_col, fc_col)
    return dict(zip(unique_df[symbol_col], unique_df[fc_col]))


def read_dmf_deg(path):
    df = pd.read_csv(path)
    require_columns(df, ["Symbol"], "DMF DEG file")
    if "log2FC" not in df.columns:
        if "log2(fc)" in df.columns:
            df["log2FC"] = df["log2(fc)"]
        else:
            raise ValueError("DMF DEG file requires either `log2FC` or `log2(fc)`.")
    return deduplicate_by_abs_fc(df, "Symbol", "log2FC")


def read_deg_csv(path, label):
    df = pd.read_csv(path)
    require_columns(df, ["Symbol", "log2FC"], label)
    return deduplicate_by_abs_fc(df, "Symbol", "log2FC")


def load_cached_deg_tables(source_results_dir):
    oxa_path = os.path.join(source_results_dir, "00_OXA_DEGs.csv")
    uru_path = os.path.join(source_results_dir, "00_URU_DEGs.csv")
    if not (os.path.exists(oxa_path) and os.path.exists(uru_path)):
        return None, None
    return (
        read_deg_csv(oxa_path, "OXA DEG table"),
        read_deg_csv(uru_path, "URU DEG table"),
    )


def extract_symbol(value):
    text = clean_symbol(value)
    if not text or text == "---":
        return ""
    if "//" in text:
        parts = [item.strip() for item in text.split("//") if item.strip()]
        if len(parts) >= 2:
            return parts[1]
        if parts:
            return parts[0]
    return text


def build_probe_to_symbol(annotation_df):
    symbol_column = None
    for column in annotation_df.columns:
        lowered = column.lower()
        if "symbol" in lowered or "gene_assignment" in lowered:
            symbol_column = column
            break
    if symbol_column is None:
        for column in annotation_df.columns:
            if "gene" in column.lower():
                symbol_column = column
                break
    if symbol_column is None:
        raise ValueError("Could not find a symbol/gene annotation column in GEO data.")

    probe_to_symbol = {}
    for _, row in annotation_df.iterrows():
        probe_id = clean_symbol(row["ID"])
        symbol = extract_symbol(row[symbol_column])
        if probe_id and symbol and symbol != "---":
            probe_to_symbol[probe_id] = symbol

    return probe_to_symbol, symbol_column


def diff_analysis(expr_df, ctrl_cols, treat_cols, probe_to_symbol, fc_cutoff, p_cutoff):
    rows = []
    log2_cutoff = np.log2(fc_cutoff)

    for probe in expr_df.index:
        ctrl_values = expr_df.loc[probe, ctrl_cols].values.astype(float)
        treat_values = expr_df.loc[probe, treat_cols].values.astype(float)
        if np.std(ctrl_values) == 0 and np.std(treat_values) == 0:
            continue

        mean_ctrl = float(np.mean(ctrl_values))
        mean_treat = float(np.mean(treat_values))
        log2_fc = mean_treat - mean_ctrl

        try:
            _, p_value = ttest_ind(treat_values, ctrl_values, equal_var=False)
        except Exception:
            continue

        if np.isnan(p_value):
            continue

        symbol = clean_symbol(probe_to_symbol.get(str(probe), ""))
        if not symbol:
            continue

        rows.append({
            "Probe": probe,
            "Symbol": symbol,
            "Ctrl_mean": mean_ctrl,
            "Treat_mean": mean_treat,
            "log2FC": log2_fc,
            "PValue": p_value,
        })

    result_df = pd.DataFrame(rows)
    if result_df.empty:
        return result_df

    result_df = result_df[
        (result_df["PValue"] < p_cutoff) &
        (result_df["log2FC"].abs() > log2_cutoff)
    ].copy()
    return deduplicate_by_abs_fc(result_df, "Symbol", "log2FC")


def load_gse(project_dir):
    try:
        import GEOparse
    except ImportError:
        raise ImportError(
            "GEOparse is required to recompute OXA/URU DEG tables. "
            "Install it or reuse cached DEG CSV files from Results/."
        )

    local_soft = os.path.join(project_dir, "GSE131963_family.soft.gz")
    if os.path.exists(local_soft):
        for keyword in ("filepath", "filename"):
            try:
                return GEOparse.get_GEO(**{keyword: local_soft, "silent": True})
            except TypeError:
                continue
            except Exception:
                pass

    return GEOparse.get_GEO(geo="GSE131963", destdir=project_dir, silent=True)


def recompute_oxa_uru_deg(project_dir, fc_cutoff, p_cutoff):
    gse = load_gse(project_dir)
    gpl = list(gse.gpls.values())[0]
    probe_to_symbol, symbol_column = build_probe_to_symbol(gpl.table)

    print("Loaded GEO annotation column: {0}".format(symbol_column))
    print("Mapped probes: {0}".format(len(probe_to_symbol)))

    expr_data = {}
    for sample_name, gsm in gse.gsms.items():
        table = gsm.table.set_index("ID_REF")
        expr_data[sample_name] = table["VALUE"].astype(float)

    expr_df = pd.DataFrame(expr_data)
    oxa_deg = diff_analysis(expr_df, CTRL_SAMPLES, OXA_SAMPLES, probe_to_symbol, fc_cutoff, p_cutoff)
    uru_deg = diff_analysis(expr_df, CTRL_SAMPLES, URU_SAMPLES, probe_to_symbol, fc_cutoff, p_cutoff)
    return oxa_deg, uru_deg


def save_deg_tables(output_dir, oxa_deg, uru_deg):
    oxa_deg.to_csv(os.path.join(output_dir, "00_OXA_DEGs.csv"), index=False)
    uru_deg.to_csv(os.path.join(output_dir, "00_URU_DEGs.csv"), index=False)


def build_gene_sets(dmf_df, oxa_deg, uru_deg):
    dmf_genes = set(normalize_symbol_series(dmf_df["Symbol"]))
    oxa_genes = set(normalize_symbol_series(oxa_deg["Symbol"]))
    uru_genes = set(normalize_symbol_series(uru_deg["Symbol"]))
    dmf_genes.discard("")
    oxa_genes.discard("")
    uru_genes.discard("")
    return dmf_genes, oxa_genes, uru_genes


def save_venn_gene_lists(output_dir, dmf_genes, oxa_genes, uru_genes):
    max_len = max(len(dmf_genes), len(oxa_genes), len(uru_genes))

    def pad(items):
        values = sorted(list(items))
        return values + [""] * (max_len - len(values))

    gene_df = pd.DataFrame({
        "DMF": pad(dmf_genes),
        "OXA": pad(oxa_genes),
        "URU": pad(uru_genes),
    })
    gene_df.to_csv(os.path.join(output_dir, "01_Venn_GeneLists.csv"), index=False)


def plot_venn_diagram(output_dir, dmf_genes, oxa_genes, uru_genes, title_suffix):
    all_three = dmf_genes & oxa_genes & uru_genes
    dmf_oxa = (dmf_genes & oxa_genes) - uru_genes
    dmf_uru = (dmf_genes & uru_genes) - oxa_genes
    oxa_uru = (oxa_genes & uru_genes) - dmf_genes
    dmf_only = dmf_genes - oxa_genes - uru_genes
    oxa_only = oxa_genes - dmf_genes - uru_genes
    uru_only = uru_genes - dmf_genes - oxa_genes

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.add_patch(Circle((-0.6, 0.4), 1.6, alpha=0.25, color=DATASET_COLORS["DMF"]))
    ax.add_patch(Circle((0.6, 0.4), 1.6, alpha=0.25, color=DATASET_COLORS["OXA"]))
    ax.add_patch(Circle((0.0, -0.5), 1.6, alpha=0.25, color=DATASET_COLORS["URU"]))

    ax.text(-1.5, 1.0, "DMF\n({0})".format(len(dmf_genes)), ha="center",
            fontsize=12, fontweight="bold", color=DATASET_COLORS["DMF"])
    ax.text(1.5, 1.0, "OXA\n({0})".format(len(oxa_genes)), ha="center",
            fontsize=12, fontweight="bold", color=DATASET_COLORS["OXA"])
    ax.text(0.0, -1.8, "URU\n({0})".format(len(uru_genes)), ha="center",
            fontsize=12, fontweight="bold", color=DATASET_COLORS["URU"])

    ax.text(-1.3, 0.0, str(len(dmf_only)), ha="center", fontsize=14, fontweight="bold")
    ax.text(1.3, 0.0, str(len(oxa_only)), ha="center", fontsize=14, fontweight="bold")
    ax.text(0.0, -1.2, str(len(uru_only)), ha="center", fontsize=14, fontweight="bold")
    ax.text(0.0, 0.9, str(len(dmf_oxa)), ha="center", fontsize=14, fontweight="bold")
    ax.text(-0.7, -0.5, str(len(dmf_uru)), ha="center", fontsize=14, fontweight="bold")
    ax.text(0.7, -0.5, str(len(oxa_uru)), ha="center", fontsize=14, fontweight="bold")
    ax.text(0.0, 0.1, str(len(all_three)), ha="center", fontsize=16,
            fontweight="bold", color="darkred")

    ax.set_xlim(-3, 3)
    ax.set_ylim(-2.8, 2.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("DEG Overlap: DMF vs OXA vs URU\n{0}".format(title_suffix),
                 fontsize=14, fontweight="bold")

    plt.savefig(os.path.join(output_dir, "01_Venn_Diagram.pdf"))
    plt.savefig(os.path.join(output_dir, "01_Venn_Diagram.png"))
    plt.close()

    return {
        "all_three": all_three,
        "dmf_oxa_only": dmf_oxa,
        "dmf_uru_only": dmf_uru,
        "oxa_uru_only": oxa_uru,
        "dmf_only": dmf_only,
        "oxa_only": oxa_only,
        "uru_only": uru_only,
    }


def build_common_fc_matrix(common_genes, dmf_df, oxa_deg, uru_deg):
    dmf_map = build_symbol_to_fc(dmf_df, "Symbol", "log2FC")
    oxa_map = build_symbol_to_fc(oxa_deg, "Symbol", "log2FC")
    uru_map = build_symbol_to_fc(uru_deg, "Symbol", "log2FC")

    rows = []
    for gene in common_genes:
        rows.append({
            "Gene": gene,
            "DMF": dmf_map.get(gene, 0.0),
            "OXA": oxa_map.get(gene, 0.0),
            "URU": uru_map.get(gene, 0.0),
        })

    if not rows:
        return pd.DataFrame(columns=["DMF", "OXA", "URU"])

    return pd.DataFrame(rows).set_index("Gene")


def plot_common_heatmap(fc_matrix, output_dir):
    fc_matrix.to_csv(os.path.join(output_dir, "02_Common_Genes_log2FC.csv"))
    if fc_matrix.empty:
        return

    ordered = fc_matrix.copy()
    ordered["mean_fc"] = ordered.mean(axis=1)
    ordered = ordered.sort_values("mean_fc", ascending=False).drop("mean_fc", axis=1)
    vmax = float(min(ordered.abs().max().max(), 6))

    if len(ordered.index) == 1:
        fig, ax = plt.subplots(figsize=(4, 3))
        sns.heatmap(ordered, cmap=HEATMAP_CMAP_R, center=0, vmin=-vmax, vmax=vmax,
                    cbar_kws={"label": "log2(FC)"}, ax=ax)
        ax.set_title("Common DEGs (n=1)", fontsize=11, fontweight="bold")
        plt.savefig(os.path.join(output_dir, "03_Common_Genes_Heatmap.pdf"))
        plt.savefig(os.path.join(output_dir, "03_Common_Genes_Heatmap.png"))
        plt.close()
        return

    grid = sns.clustermap(
        ordered,
        cmap=HEATMAP_CMAP_R,
        center=0,
        vmin=-vmax,
        vmax=vmax,
        row_cluster=True,
        col_cluster=False,
        figsize=(4.5, 12),
        linewidths=0,
        yticklabels=False,
        xticklabels=True,
        dendrogram_ratio=(0.1, 0.0),
        cbar_pos=None,
    )
    cbar_ax = grid.fig.add_axes([0.92, 0.4, 0.03, 0.2])
    colorbar = grid.fig.colorbar(grid.ax_heatmap.collections[0], cax=cbar_ax)
    colorbar.set_label("log2(FC)", fontsize=9)
    colorbar.ax.tick_params(labelsize=7)
    heatmap_position = grid.ax_heatmap.get_position()
    grid.ax_heatmap.set_position([heatmap_position.x0, heatmap_position.y0, 0.70, heatmap_position.height])
    dendrogram_position = grid.ax_row_dendrogram.get_position()
    grid.ax_row_dendrogram.set_position([
        dendrogram_position.x0,
        dendrogram_position.y0,
        dendrogram_position.width,
        heatmap_position.height,
    ])
    grid.ax_heatmap.set_title("Common DEGs (n={0})".format(len(ordered.index)),
                              fontsize=11, fontweight="bold", pad=10)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), fontsize=10, fontweight="bold")
    grid.savefig(os.path.join(output_dir, "03_Common_Genes_Heatmap.pdf"))
    grid.savefig(os.path.join(output_dir, "03_Common_Genes_Heatmap.png"))
    plt.close()


def try_import_gseapy():
    try:
        import gseapy
        return gseapy
    except ImportError:
        return None


def copy_existing_enrichment_csvs(source_results_dir, output_dir):
    copied = []
    for prefix in ("04_GO_BP", "04_GO_MF", "05_KEGG"):
        source_csv = os.path.join(source_results_dir, prefix + "_Enrichment.csv")
        target_csv = os.path.join(output_dir, prefix + "_Enrichment.csv")
        if os.path.exists(source_csv):
            shutil.copy2(source_csv, target_csv)
            copied.append(prefix)
    return copied


def clusterprofiler_dotplot(csv_path, output_prefix, title, output_dir, common_gene_count):
    enrichment_df = pd.read_csv(csv_path)
    if enrichment_df.empty:
        return

    enrichment_df = enrichment_df.sort_values("Adjusted P-value").head(15).copy()
    enrichment_df["Count"] = enrichment_df["Overlap"].map(
        lambda text: int(str(text).split("/")[0]) if "/" in str(text) else 0
    )
    enrichment_df["GeneRatio"] = enrichment_df["Count"] / max(float(common_gene_count), 1.0)
    enrichment_df["p.adjust"] = enrichment_df["Adjusted P-value"]
    enrichment_df["Description"] = enrichment_df["Term"].map(
        lambda text: text.split(" (GO:")[0] if " (GO:" in str(text) else str(text)
    )
    enrichment_df = enrichment_df.sort_values("Count", ascending=True).reset_index(drop=True)
    enrichment_df["Label"] = enrichment_df["Description"].map(
        lambda text: text[:60] + "..." if len(str(text)) > 60 else str(text)
    )

    fig_height = max(5, len(enrichment_df) * 0.45)
    fig, ax = plt.subplots(figsize=(8, fig_height))

    norm = Normalize(
        vmin=enrichment_df["p.adjust"].min(),
        vmax=enrichment_df["p.adjust"].max(),
    )
    cmap = LinearSegmentedColormap.from_list("light_rdbu", HEATMAP_COLORS)
    bubble_colors = [cmap(norm(value)) for value in enrichment_df["p.adjust"]]

    size_min = 40
    size_max = 350
    count_min = enrichment_df["Count"].min()
    count_max = enrichment_df["Count"].max()
    if count_max == count_min:
        sizes = np.array([200] * len(enrichment_df))
    else:
        sizes = size_min + (
            (enrichment_df["Count"] - count_min) /
            float(count_max - count_min)
        ) * (size_max - size_min)

    ax.scatter(
        enrichment_df["GeneRatio"],
        range(len(enrichment_df)),
        s=sizes,
        c=bubble_colors,
        edgecolors="grey",
        linewidth=0.5,
        zorder=5,
    )
    ax.set_yticks(range(len(enrichment_df)))
    ax.set_yticklabels(enrichment_df["Label"], fontsize=9)
    ax.set_xlabel("GeneRatio", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_facecolor("white")
    ax.grid(True, color="#CCCCCC", linewidth=0.5, zorder=0)

    for spine in ax.spines.values():
        spine.set_color("#666666")
        spine.set_linewidth(0.8)

    scalar = mcm.ScalarMappable(cmap=cmap, norm=norm)
    scalar.set_array([])
    cbar = plt.colorbar(scalar, ax=ax, shrink=0.35, pad=0.02, aspect=12)
    cbar.set_label("p.adjust", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    unique_counts = sorted(set(enrichment_df["Count"].tolist()))
    if len(unique_counts) >= 5:
        legend_counts = [
            unique_counts[0],
            unique_counts[len(unique_counts) // 3],
            unique_counts[(2 * len(unique_counts)) // 3],
            unique_counts[-1],
        ]
    elif len(unique_counts) >= 3:
        legend_counts = [
            unique_counts[0],
            unique_counts[len(unique_counts) // 2],
            unique_counts[-1],
        ]
    else:
        legend_counts = unique_counts
    legend_counts = sorted(set(legend_counts))

    handles = []
    for count in legend_counts:
        if count_max == count_min:
            legend_size = 200
        else:
            legend_size = size_min + ((count - count_min) / float(count_max - count_min)) * (size_max - size_min)
        handles.append(ax.scatter([], [], s=legend_size, facecolors="none",
                                  edgecolors="black", linewidth=0.8))

    if handles:
        ax.legend(
            handles,
            [str(count) for count in legend_counts],
            loc="lower right",
            fontsize=8,
            frameon=True,
            title="Count",
            title_fontsize=9,
            scatterpoints=1,
            labelspacing=1.2,
            borderpad=1,
            handletextpad=1,
            framealpha=1,
            edgecolor="#CCCCCC",
        )

    plt.savefig(os.path.join(output_dir, output_prefix + "_Dotplot.pdf"))
    plt.savefig(os.path.join(output_dir, output_prefix + "_Dotplot.png"))
    plt.close()


def run_enrichment(common_genes, source_results_dir, output_dir):
    if len(common_genes) <= 5:
        return "Skipped enrichment because there were <= 5 common genes."

    gseapy = try_import_gseapy()
    generated = []

    if gseapy is not None:
        libraries = [
            ("GO Biological Process", "GO_Biological_Process_2023", "04_GO_BP"),
            ("GO Molecular Function", "GO_Molecular_Function_2023", "04_GO_MF"),
            ("KEGG Pathway", "KEGG_2021_Human", "05_KEGG"),
        ]
        gene_list_upper = [gene.upper() for gene in common_genes]

        for title, library_name, file_prefix in libraries:
            try:
                enrichment = gseapy.enrichr(
                    gene_list=gene_list_upper,
                    gene_sets=library_name,
                    organism="mouse",
                    outdir=None,
                    no_plot=True,
                )
                result_df = enrichment.results
                result_df = result_df[result_df["Adjusted P-value"] < 0.05].head(20)
                if result_df.empty:
                    continue

                csv_path = os.path.join(output_dir, file_prefix + "_Enrichment.csv")
                result_df.to_csv(csv_path, index=False)
                clusterprofiler_dotplot(csv_path, file_prefix, title, output_dir, len(common_genes))
                generated.append(file_prefix)
            except Exception as exc:
                print("Enrichment failed for {0}: {1}".format(title, exc))

        if generated:
            return "Generated enrichment outputs for: {0}".format(", ".join(generated))

    copied = copy_existing_enrichment_csvs(source_results_dir, output_dir)
    for prefix, title in [
        ("04_GO_BP", "GO Biological Process"),
        ("04_GO_MF", "GO Molecular Function"),
        ("05_KEGG", "KEGG Pathway"),
    ]:
        csv_path = os.path.join(output_dir, prefix + "_Enrichment.csv")
        if os.path.exists(csv_path):
            clusterprofiler_dotplot(csv_path, prefix, title, output_dir, len(common_genes))

    if copied:
        return "Reused copied enrichment CSV files from Results/: {0}".format(", ".join(copied))

    return "Skipped enrichment because gseapy was unavailable and no copied CSVs were found."


def get_marker_fc(deg_df, markers_dict, symbol_col):
    result = {}
    working = deg_df.copy()
    working[symbol_col] = normalize_symbol_series(working[symbol_col])

    for th_type, markers in markers_dict.items():
        for gene in markers:
            matched = working[working[symbol_col] == gene]
            if len(matched.index) > 0:
                best_index = matched["log2FC"].abs().idxmax()
                result[gene] = {
                    "log2FC": float(matched.loc[best_index, "log2FC"]),
                    "PValue": float(matched.loc[best_index, "PValue"]) if "PValue" in matched.columns else np.nan,
                    "Th_type": th_type,
                    "is_DEG": True,
                }
            else:
                result[gene] = {
                    "log2FC": 0.0,
                    "PValue": np.nan,
                    "Th_type": th_type,
                    "is_DEG": False,
                }
    return result


def significance_stars(p_value):
    if pd.isnull(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def calc_scores(marker_fc, markers_dict):
    scores = {}
    for th_type, markers in markers_dict.items():
        fc_values = [marker_fc.get(gene, {}).get("log2FC", 0.0) for gene in markers]
        n_deg = sum(1 for gene in markers if marker_fc.get(gene, {}).get("is_DEG", False))
        scores[th_type] = {
            "mean_log2FC": float(np.mean(fc_values)),
            "n_deg": n_deg,
            "deg_ratio": n_deg / float(len(markers)),
            "n_total": len(markers),
        }
    return scores


def save_score_table(output_dir, dmf_scores, oxa_scores, uru_scores):
    rows = []
    for th_type in TH_MARKERS_MOUSE:
        rows.append({
            "Th_type": th_type,
            "DMF_meanFC": round(dmf_scores[th_type]["mean_log2FC"], 3),
            "OXA_meanFC": round(oxa_scores[th_type]["mean_log2FC"], 3),
            "URU_meanFC": round(uru_scores[th_type]["mean_log2FC"], 3),
            "DMF_nDEG": dmf_scores[th_type]["n_deg"],
            "OXA_nDEG": oxa_scores[th_type]["n_deg"],
            "URU_nDEG": uru_scores[th_type]["n_deg"],
        })

    score_df = pd.DataFrame(rows)
    score_df.to_csv(os.path.join(output_dir, "06_Th_Scores.csv"), index=False)
    return score_df


def plot_marker_heatmap(output_dir, dmf_marker_fc, oxa_marker_fc, uru_marker_fc):
    ordered_genes = []
    gene_to_type = {}
    seen = set()
    for th_type, markers in TH_MARKERS_MOUSE.items():
        for gene in markers:
            if gene not in seen:
                ordered_genes.append(gene)
                gene_to_type[gene] = th_type
                seen.add(gene)

    fc_matrix = pd.DataFrame(index=ordered_genes, columns=["DMF", "OXA", "URU"], dtype=float)
    star_matrix = pd.DataFrame(index=ordered_genes, columns=["DMF", "OXA", "URU"], dtype=object)

    for gene in ordered_genes:
        fc_matrix.loc[gene, "DMF"] = dmf_marker_fc.get(gene, {}).get("log2FC", 0.0)
        fc_matrix.loc[gene, "OXA"] = oxa_marker_fc.get(gene, {}).get("log2FC", 0.0)
        fc_matrix.loc[gene, "URU"] = uru_marker_fc.get(gene, {}).get("log2FC", 0.0)
        star_matrix.loc[gene, "DMF"] = significance_stars(dmf_marker_fc.get(gene, {}).get("PValue", np.nan))
        star_matrix.loc[gene, "OXA"] = significance_stars(oxa_marker_fc.get(gene, {}).get("PValue", np.nan))
        star_matrix.loc[gene, "URU"] = significance_stars(uru_marker_fc.get(gene, {}).get("PValue", np.nan))

    vmax = float(min(fc_matrix.abs().max().max(), 6))
    fig, axes = plt.subplots(1, 2, figsize=(10, 16),
                             gridspec_kw={"width_ratios": [3, 0.5], "wspace": 0.02})
    heatmap_ax = axes[0]
    image = heatmap_ax.imshow(fc_matrix.values, cmap=HEATMAP_CMAP_R,
                              aspect="auto", vmin=-vmax, vmax=vmax)

    for row_index in range(len(ordered_genes)):
        for col_index in range(3):
            star_label = star_matrix.iloc[row_index, col_index]
            if star_label:
                heatmap_ax.text(col_index, row_index, star_label, ha="center", va="center",
                                fontsize=7, fontweight="bold", color="black")

    heatmap_ax.set_xticks(range(3))
    heatmap_ax.set_xticklabels(["DMF", "OXA", "URU"], fontsize=10, fontweight="bold")
    heatmap_ax.set_yticks(range(len(ordered_genes)))
    heatmap_ax.set_yticklabels(ordered_genes, fontsize=7)

    th_types = [gene_to_type[gene] for gene in ordered_genes]
    previous = None
    for row_index, th_type in enumerate(th_types):
        if th_type != previous and row_index > 0:
            heatmap_ax.axhline(y=row_index - 0.5, color="white", linewidth=2)
        previous = th_type

    heatmap_ax.set_title("Th Marker log2(FC) - Mouse ACD Models\n* P<0.05, ** P<0.01, *** P<0.001",
                         fontsize=13, fontweight="bold", pad=15)

    annotation_ax = axes[1]
    annotation_ax.set_xlim(0, 1)
    annotation_ax.set_ylim(-0.5, len(ordered_genes) - 0.5)
    annotation_ax.invert_yaxis()
    for row_index, th_type in enumerate(th_types):
        annotation_ax.add_patch(plt.Rectangle(
            (0, row_index - 0.5), 1, 1,
            facecolor=TH_COLORS[th_type],
            edgecolor="white",
            linewidth=0.5,
        ))
    annotation_ax.set_xticks([])
    annotation_ax.set_yticks([])
    annotation_ax.set_xlabel("Th type", fontsize=9, fontweight="bold")

    previous = None
    starts = []
    for row_index, th_type in enumerate(th_types):
        if th_type != previous:
            starts.append((row_index, th_type))
            previous = th_type
    for start_index, th_type in starts:
        end_index = start_index
        for probe_index in range(start_index, len(th_types)):
            if th_types[probe_index] == th_type:
                end_index = probe_index
            else:
                break
        annotation_ax.text(0.5, (start_index + end_index) / 2.0, th_type,
                           ha="center", va="center", fontsize=8,
                           fontweight="bold", color="white")

    colorbar = fig.colorbar(image, ax=axes, shrink=0.3, aspect=20, pad=0.02)
    colorbar.set_label("log2(FC)")

    plt.savefig(os.path.join(output_dir, "07_Th_Marker_Heatmap.pdf"))
    plt.savefig(os.path.join(output_dir, "07_Th_Marker_Heatmap.png"))
    plt.close()


def plot_score_heatmap(output_dir, score_df):
    heatmap_df = score_df.set_index("Th_type")[["DMF_meanFC", "OXA_meanFC", "URU_meanFC"]]
    heatmap_df.columns = ["DMF", "OXA", "URU"]

    fig, ax = plt.subplots(figsize=(7, 5))
    vmax = max(abs(heatmap_df.values.min()), abs(heatmap_df.values.max()))
    image = ax.imshow(heatmap_df.values, cmap=HEATMAP_CMAP_R, aspect="auto", vmin=-vmax, vmax=vmax)

    for row_index in range(heatmap_df.shape[0]):
        for col_index in range(heatmap_df.shape[1]):
            value = heatmap_df.iloc[row_index, col_index]
            text_color = "white" if abs(value) > vmax * 0.6 else "black"
            ax.text(col_index, row_index, "{0:.2f}".format(value), ha="center", va="center",
                    fontsize=10, fontweight="bold", color=text_color)

    ax.set_xticks(range(3))
    ax.set_xticklabels(["DMF", "OXA", "URU"], fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(heatmap_df.index)))
    labels = ax.set_yticklabels(heatmap_df.index, fontsize=11, fontweight="bold")
    for label in labels:
        label.set_color(TH_COLORS.get(label.get_text(), "black"))

    ax.set_title("Th Subtype Activation Score\n(Mean log2FC of Marker Genes)",
                 fontsize=13, fontweight="bold")
    colorbar = plt.colorbar(image, ax=ax, shrink=0.8)
    colorbar.set_label("Mean log2(FC)")
    plt.savefig(os.path.join(output_dir, "08_Th_Score_Heatmap.pdf"))
    plt.savefig(os.path.join(output_dir, "08_Th_Score_Heatmap.png"))
    plt.close()


def plot_radar_overlay(output_dir, dmf_scores, oxa_scores, uru_scores):
    radar_types = ["Th1", "Th2", "Th17", "Th22", "Treg", "Tfh", "Th9"]
    angles = np.linspace(0, 2 * np.pi, len(radar_types), endpoint=False).tolist()
    angles += angles[:1]
    datasets = [
        ("DMF", dmf_scores, DATASET_COLORS["DMF"]),
        ("OXA", oxa_scores, DATASET_COLORS["OXA"]),
        ("URU", uru_scores, DATASET_COLORS["URU"]),
    ]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    for name, score_map, color in datasets:
        values = [score_map[th_type]["mean_log2FC"] for th_type in radar_types]
        values += values[:1]
        ax.plot(angles, values, "o-", linewidth=2.5, color=color, markersize=7, label=name)
        ax.fill(angles, values, alpha=0.08, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(radar_types, fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=11)
    ax.set_title("Th Subtype Comparison\nDMF vs OXA vs URU",
                 fontsize=14, fontweight="bold", pad=25)
    plt.savefig(os.path.join(output_dir, "09_Th_Radar_Overlay.pdf"))
    plt.savefig(os.path.join(output_dir, "09_Th_Radar_Overlay.png"))
    plt.close()


def plot_deg_detection_rate(output_dir, dmf_scores, oxa_scores, uru_scores):
    radar_types = ["Th1", "Th2", "Th17", "Th22", "Treg", "Tfh", "Th9"]
    datasets = [
        ("DMF", dmf_scores, DATASET_COLORS["DMF"]),
        ("OXA", oxa_scores, DATASET_COLORS["OXA"]),
        ("URU", uru_scores, DATASET_COLORS["URU"]),
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    x_axis = np.arange(len(radar_types))
    width = 0.25

    for dataset_index, (name, score_map, color) in enumerate(datasets):
        values = [score_map[th_type]["deg_ratio"] * 100.0 for th_type in radar_types]
        bars = ax.bar(x_axis + dataset_index * width, values, width, label=name,
                      color=color, edgecolor="white")
        for bar, value in zip(bars, values):
            if value > 0:
                ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() + 1,
                        "{0:.0f}%".format(value), ha="center", va="bottom",
                        fontsize=7, fontweight="bold")

    ax.set_xticks(x_axis + width)
    ax.set_xticklabels(radar_types, fontsize=11, fontweight="bold")
    ax.set_ylabel("DEG Detection Rate (%)", fontsize=11)
    ax.set_title("Th Markers Detected as DEGs", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_ylim(0, 110)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig(os.path.join(output_dir, "10_DEG_Detection_Rate.pdf"))
    plt.savefig(os.path.join(output_dir, "10_DEG_Detection_Rate.png"))
    plt.close()


def plot_polarization(output_dir, dmf_scores, oxa_scores, uru_scores):
    datasets = [
        ("DMF", dmf_scores, DATASET_COLORS["DMF"]),
        ("OXA", oxa_scores, DATASET_COLORS["OXA"]),
        ("URU", uru_scores, DATASET_COLORS["URU"]),
    ]

    fig, ax = plt.subplots(figsize=(8, 7))
    for name, score_map, color in datasets:
        th1_minus_th2 = score_map["Th1"]["mean_log2FC"] - score_map["Th2"]["mean_log2FC"]
        th17_minus_treg = score_map["Th17"]["mean_log2FC"] - score_map["Treg"]["mean_log2FC"]
        ax.scatter(th1_minus_th2, th17_minus_treg, s=300, c=color,
                   edgecolors="black", linewidth=1.5, zorder=5)
        ax.annotate(name, (th1_minus_th2, th17_minus_treg), fontsize=12,
                    fontweight="bold", xytext=(15, 10), textcoords="offset points",
                    color=color)

    ax.axhline(0, color="gray", linewidth=1, linestyle="--", alpha=0.5)
    ax.axvline(0, color="gray", linewidth=1, linestyle="--", alpha=0.5)
    xlim = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1])) * 1.3
    ylim = max(abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1])) * 1.3
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-ylim, ylim)
    ax.text(xlim * 0.7, ylim * 0.85, "Th1+Th17", ha="center", fontsize=10,
            style="italic", color="gray", alpha=0.7)
    ax.text(-xlim * 0.7, ylim * 0.85, "Th2+Th17", ha="center", fontsize=10,
            style="italic", color="gray", alpha=0.7)
    ax.text(xlim * 0.7, -ylim * 0.85, "Th1+Treg", ha="center", fontsize=10,
            style="italic", color="gray", alpha=0.7)
    ax.text(-xlim * 0.7, -ylim * 0.85, "Th2+Treg", ha="center", fontsize=10,
            style="italic", color="gray", alpha=0.7)
    ax.set_xlabel("Th1-Th2 Balance", fontsize=11, fontweight="bold")
    ax.set_ylabel("Th17-Treg Balance", fontsize=11, fontweight="bold")
    ax.set_title("Th Polarization Landscape\nMouse ACD Models",
                 fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig(os.path.join(output_dir, "11_Th_Polarization.pdf"))
    plt.savefig(os.path.join(output_dir, "11_Th_Polarization.png"))
    plt.close()


def plot_th_composition(output_dir, dmf_scores, oxa_scores, uru_scores):
    radar_types = ["Th1", "Th2", "Th17", "Th22", "Treg", "Tfh", "Th9"]
    dataset_names = ["DMF", "OXA", "URU"]
    all_scores = [dmf_scores, oxa_scores, uru_scores]

    fig, ax = plt.subplots(figsize=(8, 6))
    bottom = np.zeros(3)
    for th_type in radar_types:
        values = [max(score_map[th_type]["mean_log2FC"], 0.0) for score_map in all_scores]
        totals = [
            sum(max(score_map[item]["mean_log2FC"], 0.0) for item in radar_types)
            for score_map in all_scores
        ]
        percentages = [
            value / total * 100.0 if total > 0 else 0.0
            for value, total in zip(values, totals)
        ]
        ax.bar(dataset_names, percentages, bottom=bottom, label=th_type,
               color=TH_COLORS[th_type], edgecolor="white", linewidth=0.5)
        for index, (percentage, current_bottom) in enumerate(zip(percentages, bottom)):
            if percentage > 3:
                ax.text(index, current_bottom + percentage / 2.0, "{0:.0f}%".format(percentage),
                        ha="center", va="center", fontsize=8, fontweight="bold", color="white")
        bottom += percentages

    ax.set_ylabel("Relative Contribution (%)", fontsize=11)
    ax.set_title("Th Subtype Composition\nMouse ACD Models",
                 fontsize=13, fontweight="bold")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
    ax.set_ylim(0, 105)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.savefig(os.path.join(output_dir, "12_Th_Composition.pdf"))
    plt.savefig(os.path.join(output_dir, "12_Th_Composition.png"))
    plt.close()


def write_summary(output_dir, gene_sets, overlaps, dmf_scores, oxa_scores, uru_scores, enrichment_status):
    summary_path = os.path.join(output_dir, "Results_Summary_repaired.txt")
    with open(summary_path, "w") as handle:
        handle.write("OXAURU_FC1.5 repaired analysis summary\n")
        handle.write("=" * 60 + "\n\n")
        handle.write("DEG counts\n")
        handle.write("DMF: {0}\n".format(len(gene_sets["DMF"])))
        handle.write("OXA: {0}\n".format(len(gene_sets["OXA"])))
        handle.write("URU: {0}\n".format(len(gene_sets["URU"])))
        handle.write("Common to all three: {0}\n".format(len(overlaps["all_three"])))
        handle.write("DMF and OXA only: {0}\n".format(len(overlaps["dmf_oxa_only"])))
        handle.write("DMF and URU only: {0}\n".format(len(overlaps["dmf_uru_only"])))
        handle.write("OXA and URU only: {0}\n\n".format(len(overlaps["oxa_uru_only"])))
        handle.write("Enrichment status\n")
        handle.write(enrichment_status + "\n\n")
        handle.write("Dominant Th subtypes\n")
        for label, score_map in [("DMF", dmf_scores), ("OXA", oxa_scores), ("URU", uru_scores)]:
            ranked = sorted(score_map.items(), key=lambda item: item[1]["mean_log2FC"], reverse=True)
            top_three = ranked[:3]
            formatted = ", ".join(
                ["{0} ({1:.2f})".format(name, item["mean_log2FC"]) for name, item in top_three]
            )
            handle.write("{0}: {1}\n".format(label, formatted))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Repaired ACD comparison analysis for DMF, OXA, and URU."
    )
    parser.add_argument(
        "--project-dir",
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Project root containing DMF-FC15.csv and Results/.",
    )
    parser.add_argument(
        "--input-results-dir",
        default=None,
        help="Existing results directory to reuse copied DEG/enrichment tables. Defaults to PROJECT_DIR/Results.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for rebuilt outputs. Defaults to PROJECT_DIR/Results_rebuilt.",
    )
    parser.add_argument(
        "--fc-cutoff",
        type=float,
        default=1.5,
        help="Fold-change cutoff used for GEO-derived OXA/URU DEGs.",
    )
    parser.add_argument(
        "--p-cutoff",
        type=float,
        default=0.05,
        help="P-value cutoff used for GEO-derived OXA/URU DEGs.",
    )
    parser.add_argument(
        "--recompute-degs",
        action="store_true",
        help="Ignore copied DEG tables and recompute OXA/URU DEGs from GEO data.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    project_dir = os.path.abspath(args.project_dir)
    source_results_dir = os.path.abspath(
        args.input_results_dir or os.path.join(project_dir, "Results")
    )
    output_dir = os.path.abspath(
        args.output_dir or os.path.join(project_dir, "Results_rebuilt")
    )

    ensure_dir(output_dir)

    dmf_path = os.path.join(project_dir, "DMF-FC15.csv")
    dmf_df = read_dmf_deg(dmf_path)

    cached_oxa, cached_uru = load_cached_deg_tables(source_results_dir)
    if not args.recompute_degs and cached_oxa is not None and cached_uru is not None:
        print("Using copied DEG tables from {0}".format(source_results_dir))
        oxa_deg = cached_oxa
        uru_deg = cached_uru
    else:
        print("Recomputing OXA/URU DEG tables from GEO data")
        oxa_deg, uru_deg = recompute_oxa_uru_deg(project_dir, args.fc_cutoff, args.p_cutoff)

    save_deg_tables(output_dir, oxa_deg, uru_deg)

    dmf_genes, oxa_genes, uru_genes = build_gene_sets(dmf_df, oxa_deg, uru_deg)
    gene_sets = {"DMF": dmf_genes, "OXA": oxa_genes, "URU": uru_genes}
    save_venn_gene_lists(output_dir, dmf_genes, oxa_genes, uru_genes)
    overlaps = plot_venn_diagram(
        output_dir,
        dmf_genes,
        oxa_genes,
        uru_genes,
        "(FC>{0}, P<{1})".format(args.fc_cutoff, args.p_cutoff),
    )

    common_genes = sorted(list(overlaps["all_three"]))
    common_fc_matrix = build_common_fc_matrix(common_genes, dmf_df, oxa_deg, uru_deg)
    plot_common_heatmap(common_fc_matrix, output_dir)

    enrichment_status = run_enrichment(common_genes, source_results_dir, output_dir)
    print(enrichment_status)

    dmf_marker_fc = get_marker_fc(dmf_df, TH_MARKERS_MOUSE, "Symbol")
    oxa_marker_fc = get_marker_fc(oxa_deg, TH_MARKERS_MOUSE, "Symbol")
    uru_marker_fc = get_marker_fc(uru_deg, TH_MARKERS_MOUSE, "Symbol")

    dmf_scores = calc_scores(dmf_marker_fc, TH_MARKERS_MOUSE)
    oxa_scores = calc_scores(oxa_marker_fc, TH_MARKERS_MOUSE)
    uru_scores = calc_scores(uru_marker_fc, TH_MARKERS_MOUSE)

    score_df = save_score_table(output_dir, dmf_scores, oxa_scores, uru_scores)
    plot_marker_heatmap(output_dir, dmf_marker_fc, oxa_marker_fc, uru_marker_fc)
    plot_score_heatmap(output_dir, score_df)
    plot_radar_overlay(output_dir, dmf_scores, oxa_scores, uru_scores)
    plot_deg_detection_rate(output_dir, dmf_scores, oxa_scores, uru_scores)
    plot_polarization(output_dir, dmf_scores, oxa_scores, uru_scores)
    plot_th_composition(output_dir, dmf_scores, oxa_scores, uru_scores)
    write_summary(output_dir, gene_sets, overlaps, dmf_scores, oxa_scores, uru_scores, enrichment_status)

    print("=" * 60)
    print("Repaired analysis completed")
    print("Output directory: {0}".format(output_dir))
    print("DMF DEGs: {0}".format(len(dmf_genes)))
    print("OXA DEGs: {0}".format(len(oxa_genes)))
    print("URU DEGs: {0}".format(len(uru_genes)))
    print("Common DEGs: {0}".format(len(common_genes)))


if __name__ == "__main__":
    main()
