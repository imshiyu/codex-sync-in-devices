from __future__ import annotations

from pathlib import Path
from textwrap import shorten

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


matplotlib.use("Agg")

INPUT = Path(r"E:\project\黑质单细胞\星胶差异表格.tsv")
OUT = Path("p_cresol_snRNA_analysis") / "astrocyte_deg_table_enrichment"

ASTRO_CLUSTER = "Astrocytes_sp_PD_vs_Astrocytes_sp_Ctrl"
LIBRARIES = [
    "GO_Biological_Process_2023",
    "KEGG_2019_Mouse",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]

FOCUS_PATTERNS = [
    "aging",
    "ageing",
    "senescence",
    "senescent",
    "parkinson",
    "neurodegeneration",
    "alzheimer",
    "huntington",
    "amyotrophic",
    "prion",
    "dopaminergic",
    "synapse",
    "synaptic",
    "mitochond",
    "oxidative phosphorylation",
    "electron transport",
    "respiratory chain",
    "autophagy",
    "mitophagy",
    "lysosome",
    "phagosome",
    "proteasome",
    "apoptosis",
    "ferroptosis",
    "necroptosis",
    "inflamm",
    "immune",
    "interferon",
    "interleukin",
    "cytokine",
    "chemokine",
    "nf-kappa",
    "tnf",
    "mapk",
    "pi3k",
    "akt",
    "mtor",
    "foxo",
    "calcium",
    "trp",
]


def import_gseapy():
    try:
        import gseapy as gp
    except ImportError as exc:
        raise RuntimeError("gseapy is required for Enrichr enrichment.") from exc
    return gp


def safe_name(text: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in str(text)).strip("_")


def clean_term(term: str) -> str:
    return str(term).replace("_", " ").rsplit(" (GO:", 1)[0].strip()


def p_to_score(pvalue: pd.Series) -> pd.Series:
    return -np.log10(pd.to_numeric(pvalue, errors="coerce").clip(lower=1e-300))


def load_astrocyte_degs() -> pd.DataFrame:
    df = pd.read_csv(INPUT, sep="\t")
    required = {"gene", "avg_log2FC", "p_val", "p_val_adj", "cluster", "type"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    astro = df[df["cluster"].astype(str) == ASTRO_CLUSTER].copy()
    if astro.empty:
        raise ValueError(f"No rows found for cluster {ASTRO_CLUSTER}")
    astro["gene"] = astro["gene"].astype(str)
    astro["avg_log2FC"] = pd.to_numeric(astro["avg_log2FC"], errors="coerce")
    astro["p_val"] = pd.to_numeric(astro["p_val"], errors="coerce")
    astro["p_val_adj"] = pd.to_numeric(astro["p_val_adj"], errors="coerce")
    astro = astro.dropna(subset=["gene", "avg_log2FC", "p_val"])
    return astro


def gene_sets_from_type(astro: pd.DataFrame) -> dict[str, list[str]]:
    gene_sets = {
        "Astrocyte_up_in_PD": astro.loc[astro["type"].astype(str) == "up", "gene"].tolist(),
        "Astrocyte_down_in_PD": astro.loc[astro["type"].astype(str) == "down", "gene"].tolist(),
    }
    # Fallback for tables without an explicit type label.
    if len(gene_sets["Astrocyte_up_in_PD"]) < 10:
        gene_sets["Astrocyte_up_in_PD"] = astro.query("avg_log2FC > 0 and p_val_adj < 0.05")["gene"].tolist()
    if len(gene_sets["Astrocyte_down_in_PD"]) < 10:
        gene_sets["Astrocyte_down_in_PD"] = astro.query("avg_log2FC < 0 and p_val_adj < 0.05")["gene"].tolist()

    # Last fallback: use the most shifted genes by nominal P value and fold change.
    if len(gene_sets["Astrocyte_up_in_PD"]) < 10:
        gene_sets["Astrocyte_up_in_PD"] = (
            astro[astro["avg_log2FC"] > 0]
            .sort_values(["p_val", "avg_log2FC"], ascending=[True, False])
            .head(250)["gene"]
            .tolist()
        )
    if len(gene_sets["Astrocyte_down_in_PD"]) < 10:
        gene_sets["Astrocyte_down_in_PD"] = (
            astro[astro["avg_log2FC"] < 0]
            .assign(abs_lfc=lambda x: x["avg_log2FC"].abs())
            .sort_values(["p_val", "abs_lfc"], ascending=[True, False])
            .head(250)["gene"]
            .tolist()
        )
    return {key: sorted(set(value)) for key, value in gene_sets.items()}


def save_gene_lists(astro: pd.DataFrame, gene_sets: dict[str, list[str]]) -> None:
    astro.to_csv(OUT / "astrocyte_DEG_table_used.csv", index=False, encoding="utf-8-sig")
    summary_rows = []
    for label, genes in gene_sets.items():
        direction = "up" if "up" in label else "down"
        pd.Series(genes, name="gene").to_csv(OUT / f"{label}_genes.txt", index=False, header=False)
        sub = astro[astro["gene"].isin(genes)].copy()
        sub.to_csv(OUT / f"{label}_DEG_table.csv", index=False, encoding="utf-8-sig")
        summary_rows.append({"label": label, "direction": direction, "n_genes": len(genes)})
    pd.DataFrame(summary_rows).to_csv(OUT / "astrocyte_DEG_gene_list_summary.csv", index=False, encoding="utf-8-sig")


def run_enrichr_for_sets(gene_sets: dict[str, list[str]], background: list[str]) -> pd.DataFrame:
    gp = import_gseapy()
    all_results = []
    log_rows = []
    for label, genes in gene_sets.items():
        label_dir = OUT / safe_name(label)
        label_dir.mkdir(parents=True, exist_ok=True)
        for library in LIBRARIES:
            try:
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=library,
                    organism="mouse",
                    outdir=str(label_dir / safe_name(library)),
                    cutoff=1.0,
                    no_plot=True,
                    background=background,
                )
                if enr is None or enr.results is None or enr.results.empty:
                    log_rows.append({"label": label, "library": library, "status": "empty", "message": ""})
                    continue
                res = enr.results.copy()
                res.insert(0, "label", label)
                res.insert(1, "library", library)
                res["direction"] = "Up in PD" if "up" in label else "Down in PD"
                all_results.append(res)
                res.to_csv(OUT / f"enrichment_{safe_name(label)}_{safe_name(library)}.csv", index=False, encoding="utf-8-sig")
                log_rows.append({"label": label, "library": library, "status": "ok", "message": f"{len(res)} rows"})
            except Exception as exc:
                message = str(exc)
                (label_dir / f"{safe_name(library)}_ERROR.txt").write_text(message, encoding="utf-8")
                log_rows.append({"label": label, "library": library, "status": "error", "message": message})
    pd.DataFrame(log_rows).to_csv(OUT / "enrichr_run_log.csv", index=False, encoding="utf-8-sig")
    if not all_results:
        return pd.DataFrame()
    combined = pd.concat(all_results, ignore_index=True)
    combined["Term_clean"] = combined["Term"].map(clean_term)
    combined["minus_log10_p"] = p_to_score(combined["P-value"])
    combined["minus_log10_fdr"] = p_to_score(combined["Adjusted P-value"])
    combined.to_csv(OUT / "astrocyte_enrichment_all_libraries.csv", index=False, encoding="utf-8-sig")
    return combined


def top_terms(df: pd.DataFrame, label: str, library: str, top_n: int = 12) -> pd.DataFrame:
    sub = df[(df["label"] == label) & (df["library"] == library)].copy()
    if sub.empty:
        return sub
    sub["Adjusted P-value"] = pd.to_numeric(sub["Adjusted P-value"], errors="coerce")
    sub["P-value"] = pd.to_numeric(sub["P-value"], errors="coerce")
    sub = sub.dropna(subset=["P-value"])
    return sub.sort_values(["Adjusted P-value", "P-value"]).head(top_n).sort_values("minus_log10_p")


def plot_bar(plot_df: pd.DataFrame, title: str, filename: str, color: str) -> None:
    if plot_df.empty:
        return
    y = np.arange(len(plot_df))
    fig, ax = plt.subplots(figsize=(9.2, max(4.6, 0.42 * len(plot_df) + 1.2)))
    ax.barh(y, plot_df["minus_log10_p"], color=color, edgecolor="black", linewidth=0.35, height=0.72)
    labels = [shorten(term, width=58, placeholder="...") for term in plot_df["Term_clean"]]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    xmax = max(2.0, float(plot_df["minus_log10_p"].max()))
    ax.set_xlim(0, np.ceil(xmax) * 1.15)
    for yi, value in zip(y, plot_df["minus_log10_p"]):
        ax.text(value + xmax * 0.015, yi, f"{value:.1f}", va="center", fontsize=9.5)
    ax.set_xlabel("-log10(P value)", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=13.5, fontweight="bold")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / f"{filename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_all_top_bars(enrichment: pd.DataFrame) -> None:
    top_rows = []
    for label, color in [("Astrocyte_up_in_PD", "#c93c3c"), ("Astrocyte_down_in_PD", "#386cb0")]:
        direction = "Up in PD" if "up" in label else "Down in PD"
        for library in LIBRARIES:
            top = top_terms(enrichment, label, library, top_n=12)
            if top.empty:
                continue
            top_rows.append(top.assign(plot_group=f"{label}_{library}"))
            plot_bar(
                top,
                f"{direction}: {library} enrichment",
                f"{safe_name(label)}_{safe_name(library)}_top{len(top)}",
                color,
            )
    if top_rows:
        pd.concat(top_rows, ignore_index=True).to_csv(
            OUT / "astrocyte_enrichment_top_terms_by_library.csv",
            index=False,
            encoding="utf-8-sig",
        )


def plot_focus_bars(enrichment: pd.DataFrame) -> None:
    pattern = "|".join(FOCUS_PATTERNS)
    focus = enrichment[enrichment["Term"].str.contains(pattern, case=False, na=False, regex=True)].copy()
    focus.to_csv(OUT / "astrocyte_PD_aging_focused_enrichment_all.csv", index=False, encoding="utf-8-sig")
    focus_rows = []
    for label, color in [("Astrocyte_up_in_PD", "#c93c3c"), ("Astrocyte_down_in_PD", "#386cb0")]:
        direction = "Up in PD" if "up" in label else "Down in PD"
        top = (
            focus[focus["label"] == label]
            .sort_values(["Adjusted P-value", "P-value"])
            .head(14)
            .sort_values("minus_log10_p")
        )
        if top.empty:
            continue
        focus_rows.append(top.assign(plot_group=label))
        plot_bar(
            top,
            f"Astrocyte {direction}: PD/aging-focused enrichment",
            f"{safe_name(label)}_PD_aging_focused_top{len(top)}",
            color,
        )
    if focus_rows:
        pd.concat(focus_rows, ignore_index=True).to_csv(
            OUT / "astrocyte_PD_aging_focused_enrichment_top_terms.csv",
            index=False,
            encoding="utf-8-sig",
        )


def write_readme(enrichment: pd.DataFrame, gene_sets: dict[str, list[str]]) -> None:
    lines = [
        "# Astrocyte DEG table enrichment",
        "",
        f"Input: `{INPUT}`",
        f"Cluster used: `{ASTRO_CLUSTER}`",
        "",
        "Gene sets:",
    ]
    for label, genes in gene_sets.items():
        lines.append(f"- `{label}`: {len(genes)} genes")
    lines.extend(
        [
            "",
            "Method:",
            "- Enrichr over-representation analysis with Mouse libraries.",
            f"- Libraries: {', '.join(LIBRARIES)}.",
            "- Up and down genes are analyzed separately.",
            "- Bar plots use `-log10(P value)`; adjusted P values are retained in CSV outputs.",
            "",
            "Interpretation:",
            "- `-log10(P value) > 1.30` corresponds to nominal P < 0.05.",
            "- For stronger inference, prioritize rows with `Adjusted P-value < 0.05`.",
            "- snRNA cell-level differential results should be treated as descriptive/exploratory without biological replicates.",
            "",
            f"Total enrichment rows: {len(enrichment)}",
        ]
    )
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    astro = load_astrocyte_degs()
    gene_sets = gene_sets_from_type(astro)
    background = sorted(astro["gene"].dropna().astype(str).unique())
    save_gene_lists(astro, gene_sets)
    enrichment = run_enrichr_for_sets(gene_sets, background)
    if enrichment.empty:
        write_readme(enrichment, gene_sets)
        raise RuntimeError("No enrichment results were returned. Check enrichr_run_log.csv for details.")
    plot_all_top_bars(enrichment)
    plot_focus_bars(enrichment)
    write_readme(enrichment, gene_sets)


if __name__ == "__main__":
    main()
