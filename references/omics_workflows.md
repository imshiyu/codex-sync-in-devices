# Omics Workflow Memory

Use these templates when Shiyu asks to rerun transcriptomics, snRNA-seq, non-targeted metabolomics, KEGG/GO enrichment, volcano plots, heatmaps, dotplots, UMAP/tSNE, or figures in the same style as previous p-cresol / PD work.

## Main Entry

List bundled tasks:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py --list
```

Run serum non-targeted metabolomics:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py serum-metabolomics --source E:\project\血清非靶 --out E:\project\血清非靶_Ctrl_vs_PD_reanalysis
```

Run midbrain non-targeted metabolomics:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py midbrain-metabolomics --source E:\project\中脑非靶 --out E:\project\中脑非靶_Ctrl_vs_PD_reanalysis_drop_Ctrl2_Ctrl3_Ctrl8_PD3_PD4_PD5
```

Scaffold snRNA-seq templates for a new project:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py snrna-template --out E:\project\new_snrna_analysis
```

Run black-substantia bulk RNA-seq FASTQ audit and same-style PPT workflow:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py midbrain-bulk-rnaseq --out H:\全转录组\黑质_RNAseq_reanalysis
```

Run the repaired skin ACD RNA-seq comparison template:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py skin-acd-rnaseq --project-dir <project-root> --output-dir <output-dir>
```

The skin ACD repaired template expects its historical project inputs, especially `DMF-FC15.csv` and either `Results/00_OXA_DEGs.csv` + `Results/00_URU_DEGs.csv` or GEO-derived files that allow recomputation. The current `H:\全转录组\皮肤有参RNA-seq` directory contains FASTQ files only, so use this template after adding the required expression/DEG inputs.

When a count/TPM matrix is available, pass it through to the template script:

```powershell
D:\Anaconda\python.exe scripts\run_omics_analysis.py midbrain-bulk-rnaseq --out H:\全转录组\黑质_RNAseq_reanalysis --matrix H:\path\to\count_matrix.csv
```

## Historical Style

Preserve the prior visual conventions unless Shiyu asks to redesign:

- Save both PNG and PDF for final figures when the template already does so.
- Use white backgrounds, compact scientific layouts, high DPI, and direct figure filenames such as `PCA score plot.png`, `volcano.png`, `heatmap.png`, `KEGG Enrichment bubble.png`, `main_umap.png`, `major_marker_dotplot.png`.
- For metabolomics, keep the Biotree-like output folder structure: `02.Statistical Analysis`, `03.Hierarchical Clustering Analysis`, `04.Boxplot Analysis`, `07.Correlation Analysis`, `10.Enrichment Analysis`, `14.ROC Curve`.
- For snRNA-seq, keep `figures/`, `tables/`, `h5ad/`, `enrichment/`, and Chinese summary reports.
- Keep the project-specific caution that cell-level snRNA tests are descriptive when Ctrl and PD each have one library.

## Bundled Templates

- `references/analysis_templates/metabolomics/rerun_serum_nontarget_ctrl_pd.py`
- `references/analysis_templates/metabolomics/rerun_midbrain_nontarget_ctrl_pd.py`
- `references/analysis_templates/metabolomics/run_filtered_nontarget_enrichment.py`
- `references/analysis_templates/transcriptomics/p_cresol_snrna_full_analysis.py`
- `references/analysis_templates/transcriptomics/analyze_midbrain_bulk_rnaseq.py`
- `references/analysis_templates/transcriptomics/skin_acd_analysis_repaired.py`
- `references/analysis_templates/transcriptomics/plot_astrocyte_b2m_heatmap_violin.py`
- `references/analysis_templates/transcriptomics/plot_astrocyte_refined_marker_dotplot.py`
- `references/analysis_templates/transcriptomics/plot_astrocyte_refined_top_marker_dotplot.py`
- `references/analysis_templates/transcriptomics/run_astrocyte_deg_table_enrichment.py`

## Dependency Hints

Metabolomics templates expect Python packages commonly used in the prior environment: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `scikit-learn`, and `openpyxl`.

snRNA-seq templates additionally expect `scanpy`, `anndata`, `gseapy`, `statsmodels`, and related scientific Python packages.
