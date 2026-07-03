from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_astrocyte_deg_table_enrichment import (
    ASTRO_CLUSTER,
    FOCUS_PATTERNS,
    INPUT,
    LIBRARIES,
    clean_term,
    import_gseapy,
    plot_bar,
    p_to_score,
    safe_name,
)


matplotlib.use("Agg")

OUT = Path("p_cresol_snRNA_analysis") / "astrocyte_nominal_p_enrichment"
P_CUTOFF = 0.05

CGAS_STING_PATTERNS = [
    "cgas",
    "sting",
    "tmem173",
    "mb21d1",
    "cytosolic dna",
    "cytosolic_dna",
    "dna sensing",
    "response to dna",
    "innate immune",
    "interferon",
    "type i interferon",
    "ddx58",
    "ifih1",
]


def load_nominal_sets() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    df = pd.read_csv(INPUT, sep="\t")
    astro = df[df["cluster"].astype(str) == ASTRO_CLUSTER].copy()
    astro["avg_log2FC"] = pd.to_numeric(astro["avg_log2FC"], errors="coerce")
    astro["p_val"] = pd.to_numeric(astro["p_val"], errors="coerce")
    astro = astro.dropna(subset=["gene", "avg_log2FC", "p_val"])
    up = astro[(astro["avg_log2FC"] > 0) & (astro["p_val"] < P_CUTOFF)]["gene"].astype(str).tolist()
    down = astro[(astro["avg_log2FC"] < 0) & (astro["p_val"] < P_CUTOFF)]["gene"].astype(str).tolist()
    return astro, {
        "Astrocyte_nominal_p_up_in_PD": sorted(set(up)),
        "Astrocyte_nominal_p_down_in_PD": sorted(set(down)),
    }


def run_enrichment(gene_sets: dict[str, list[str]], background: list[str]) -> pd.DataFrame:
    gp = import_gseapy()
    rows = []
    logs = []
    for label, genes in gene_sets.items():
        pd.Series(genes, name="gene").to_csv(OUT / f"{label}_genes.txt", index=False, header=False)
        for library in LIBRARIES:
            try:
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=library,
                    organism="mouse",
                    outdir=str(OUT / safe_name(label) / safe_name(library)),
                    cutoff=1.0,
                    no_plot=True,
                    background=background,
                )
                if enr is not None and enr.results is not None and not enr.results.empty:
                    res = enr.results.copy()
                    res.insert(0, "label", label)
                    res.insert(1, "library", library)
                    rows.append(res)
                    logs.append({"label": label, "library": library, "status": "ok", "message": len(res)})
                else:
                    logs.append({"label": label, "library": library, "status": "empty", "message": ""})
            except Exception as exc:
                logs.append({"label": label, "library": library, "status": "error", "message": str(exc)})
    pd.DataFrame(logs).to_csv(OUT / "enrichr_run_log.csv", index=False, encoding="utf-8-sig")
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    out["Term_clean"] = out["Term"].map(clean_term)
    out["minus_log10_p"] = p_to_score(out["P-value"])
    out["minus_log10_fdr"] = p_to_score(out["Adjusted P-value"])
    out.to_csv(OUT / "astrocyte_nominal_p_enrichment_all_libraries.csv", index=False, encoding="utf-8-sig")
    return out


def save_target_matches(enrichment: pd.DataFrame) -> None:
    pattern = "|".join(CGAS_STING_PATTERNS)
    matches = enrichment[
        enrichment["Term"].str.contains(pattern, case=False, regex=True, na=False)
        | enrichment["Genes"].str.contains(pattern, case=False, regex=True, na=False)
    ].copy()
    matches = matches.sort_values(["label", "Adjusted P-value", "P-value"])
    matches.to_csv(OUT / "cGAS_STING_interferon_related_enrichment_matches.csv", index=False, encoding="utf-8-sig")

    for label, color in [
        ("Astrocyte_nominal_p_up_in_PD", "#c93c3c"),
        ("Astrocyte_nominal_p_down_in_PD", "#386cb0"),
    ]:
        top = matches[matches["label"] == label].sort_values(["P-value"]).head(12).sort_values("minus_log10_p")
        direction = "Up in PD" if "up" in label else "Down in PD"
        plot_bar(
            top,
            f"Nominal-P astrocyte {direction}: cGAS/STING/interferon-related terms",
            f"{label}_cGAS_STING_interferon_related_top{len(top)}",
            color,
        )


def save_focus_top(enrichment: pd.DataFrame) -> None:
    pattern = "|".join(FOCUS_PATTERNS)
    focus = enrichment[enrichment["Term"].str.contains(pattern, case=False, regex=True, na=False)].copy()
    focus.to_csv(OUT / "astrocyte_nominal_p_PD_aging_focused_all.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    astro, gene_sets = load_nominal_sets()
    background = sorted(astro["gene"].astype(str).unique())
    pd.DataFrame(
        [{"label": label, "n_genes": len(genes)} for label, genes in gene_sets.items()]
    ).to_csv(OUT / "nominal_p_gene_list_summary.csv", index=False, encoding="utf-8-sig")
    enrichment = run_enrichment(gene_sets, background)
    if enrichment.empty:
        raise RuntimeError("No enrichment results returned.")
    save_target_matches(enrichment)
    save_focus_top(enrichment)
    (OUT / "README.md").write_text(
        "\n".join(
            [
                "# Astrocyte nominal-P enrichment",
                "",
                f"Input: `{INPUT}`",
                f"Cluster: `{ASTRO_CLUSTER}`",
                f"Nominal cutoff: `p_val < {P_CUTOFF}` with direction from `avg_log2FC`.",
                "",
                "This is exploratory and not FDR-controlled at the DEG-selection step.",
                "Use it only to inspect whether weak pathway-level signals appear under a looser gene filter.",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
