from __future__ import annotations

import math
import sys
import traceback
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(r"E:\project\血清非靶_Ctrl_vs_PD_reanalysis")
ENRICH_DIR = (
    BASE_DIR
    / "10.Enrichment Analysis"
    / "Ctrl vs PD"
    / "extended_enrichment"
)
OUT_DIR = BASE_DIR / "10.Enrichment Analysis" / "PD_vs_Ctrl_enrichment_count_legend"


def wrap_label(text: str, width: int = 42) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False))


def bubble_size(count: float) -> float:
    return 70 + 70 * float(count)


def legend_counts(values: pd.Series) -> list[int]:
    counts = sorted({int(v) for v in values.dropna().tolist()})
    if len(counts) <= 4:
        return counts
    mid = counts[len(counts) // 2]
    return sorted({counts[0], mid, counts[-1]})


def plot_kegg_bubble(
    df: pd.DataFrame,
    query_group: str,
    title: str,
    filename_stem: str,
    top_n: int = 15,
) -> None:
    sub = df.loc[df["QueryGroup"].eq(query_group)].copy()
    if sub.empty:
        return

    sub["Adjusted P-value"] = pd.to_numeric(sub["Adjusted P-value"], errors="coerce")
    sub["P-value"] = pd.to_numeric(sub["P-value"], errors="coerce")
    sub["Rich_factor"] = pd.to_numeric(sub["Rich_factor"], errors="coerce")
    sub["Overlap_count"] = pd.to_numeric(sub["Overlap_count"], errors="coerce")
    sub = sub.dropna(subset=["Description", "Rich_factor", "Overlap_count", "P-value"])
    sub["plot_p"] = sub["Adjusted P-value"].where(sub["Adjusted P-value"].gt(0), sub["P-value"])
    sub["plot_p"] = sub["plot_p"].clip(lower=1e-300)
    sub["minus_log10_adj_p"] = -sub["plot_p"].map(math.log10)

    sub = sub.sort_values(["Adjusted P-value", "P-value", "Overlap_count"], ascending=[True, True, False])
    sub = sub.head(top_n).copy()
    sub = sub.iloc[::-1].reset_index(drop=True)
    sub["Description_wrapped"] = sub["Description"].map(wrap_label)

    height = max(4.8, 0.42 * len(sub) + 1.7)
    fig, ax = plt.subplots(figsize=(8.8, height), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    scatter = ax.scatter(
        sub["Rich_factor"],
        sub["Description_wrapped"],
        s=sub["Overlap_count"].map(bubble_size),
        c=sub["minus_log10_adj_p"],
        cmap="viridis_r",
        edgecolor="#333333",
        linewidth=0.55,
        alpha=0.92,
    )

    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xlabel("Rich factor", fontsize=11)
    ax.set_ylabel("")
    ax.grid(axis="x", color="#dddddd", linewidth=0.7)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.tick_params(axis="x", labelsize=9)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")

    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.045)
    colorbar.set_label("-log10(adjusted P)", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)

    handles = [
        ax.scatter(
            [],
            [],
            s=bubble_size(count),
            facecolor="#bdbdbd",
            edgecolor="#333333",
            linewidth=0.55,
            alpha=0.92,
            label=str(count),
        )
        for count in legend_counts(sub["Overlap_count"])
    ]
    size_legend = ax.legend(
        handles=handles,
        title="Overlap count",
        frameon=True,
        loc="lower right",
        bbox_to_anchor=(1.0, 0.02),
        fontsize=8,
        title_fontsize=9,
        borderpad=0.7,
        labelspacing=0.9,
    )
    size_legend.get_frame().set_facecolor("white")
    size_legend.get_frame().set_edgecolor("#cccccc")
    ax.add_artist(size_legend)

    fig.tight_layout()
    png = OUT_DIR / f"{filename_stem}.png"
    pdf = OUT_DIR / f"{filename_stem}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def plot_da_score(df: pd.DataFrame) -> None:
    sub = df.loc[df["QueryGroup"].eq("All_DEGs")].copy()
    if sub.empty:
        return

    for col in ["Adjusted P-value", "P-value", "DA_score", "Overlap_count"]:
        sub[col] = pd.to_numeric(sub[col], errors="coerce")
    sub = sub.dropna(subset=["Description", "DA_score", "P-value", "Overlap_count"])
    sub["plot_p"] = sub["Adjusted P-value"].where(sub["Adjusted P-value"].gt(0), sub["P-value"])
    sub["plot_p"] = sub["plot_p"].clip(lower=1e-300)
    sub["minus_log10_adj_p"] = -sub["plot_p"].map(math.log10)
    sub["abs_da"] = sub["DA_score"].abs()
    sub = sub.sort_values(["abs_da", "Adjusted P-value", "Overlap_count"], ascending=[False, True, False])
    sub = sub.head(15).copy()
    sub = sub.sort_values("DA_score").reset_index(drop=True)
    sub["Description_wrapped"] = sub["Description"].map(wrap_label)

    fig, ax = plt.subplots(figsize=(9.0, max(4.8, 0.42 * len(sub) + 1.6)), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    colors = ["#3b6fb6" if value < 0 else "#c44e52" for value in sub["DA_score"]]
    ax.barh(sub["Description_wrapped"], sub["DA_score"], color=colors, alpha=0.9)
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_title("KEGG directional activity score (PD vs Ctrl)", fontsize=13, pad=12)
    ax.set_xlabel("DA score  (positive = PD-up; negative = PD-down)", fontsize=10.5)
    ax.set_ylabel("")
    ax.grid(axis="x", color="#dddddd", linewidth=0.7)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.tick_params(axis="x", labelsize=9)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "PD_vs_Ctrl_KEGG_DA_score.png", bbox_inches="tight")
    fig.savefig(OUT_DIR / "PD_vs_Ctrl_KEGG_DA_score.pdf", bbox_inches="tight")
    plt.close(fig)


def write_readme(kegg: pd.DataFrame) -> None:
    summary = (
        kegg.groupby("QueryGroup", dropna=False)
        .agg(
            pathways=("Description", "count"),
            fdr_0_05=("Adjusted P-value", lambda x: int((pd.to_numeric(x, errors="coerce") <= 0.05).sum())),
        )
        .reset_index()
    )
    summary_text = summary.to_string(index=False)

    readme = f"""# PD vs Ctrl serum KEGG enrichment plots with count legend

Input table:

`{ENRICH_DIR / "KEGG_pathway_enrichment_extended.csv"}`

This redraw keeps the existing enrichment statistics and fixes the bubble-size legend.

- X axis: Rich factor
- Dot color: -log10(adjusted P)
- Dot size: Overlap count
- Direction: Up means higher in PD than Ctrl; Down means lower in PD than Ctrl.

Summary:

{summary_text}

Outputs:

- `PD_vs_Ctrl_KEGG_All_DEGs_bubble_count_legend.png/pdf`
- `PD_vs_Ctrl_KEGG_Up_bubble_count_legend.png/pdf`
- `PD_vs_Ctrl_KEGG_Down_bubble_count_legend.png/pdf`
- `PD_vs_Ctrl_KEGG_DA_score.png/pdf`
"""
    (OUT_DIR / "README_PD_vs_Ctrl_enrichment_count_legend.md").write_text(readme, encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    kegg_path = ENRICH_DIR / "KEGG_pathway_enrichment_extended.csv"
    kegg = pd.read_csv(kegg_path)

    plot_kegg_bubble(
        kegg,
        query_group="All_DEGs",
        title="KEGG enrichment of differential serum metabolites (PD vs Ctrl)",
        filename_stem="PD_vs_Ctrl_KEGG_All_DEGs_bubble_count_legend",
    )
    plot_kegg_bubble(
        kegg,
        query_group="Up",
        title="KEGG enrichment of PD-up serum metabolites",
        filename_stem="PD_vs_Ctrl_KEGG_Up_bubble_count_legend",
    )
    plot_kegg_bubble(
        kegg,
        query_group="Down",
        title="KEGG enrichment of PD-down serum metabolites",
        filename_stem="PD_vs_Ctrl_KEGG_Down_bubble_count_legend",
    )
    plot_da_score(kegg)
    write_readme(kegg)
    print(f"Wrote updated enrichment plots to: {OUT_DIR}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc(file=sys.stderr)
        raise
