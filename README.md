# Age-dependent FOXO3 rewiring of B-cell immunity shapes outcomes in viral acute respiratory distress syndrome

**Author:** Licheng Song  
**Institution:** Chinese PLA General Hospital (Eighth Medical Center) · Capital Institute of Pediatrics  
**Prior study:** [Lung-Repair-Niches-in-Pediatric-ARDS](https://github.com/songlicheng666/Lung-Repair-Niches-in-Pediatric-ARDS)

---

## Overview

This repository contains all analysis scripts for the manuscript investigating age-dependent FOXO3-driven B-cell immune rewiring in virus-associated ARDS. Three complementary multi-omics platforms are integrated:

1. **Plasma proteomics** (DIA-MS, n=67 samples) — FOXO3-centered prognostic protein panel (AUC=0.859)
2. **Peripheral blood scRNA-seq + BCR V(D)J sequencing** (10x Genomics 5′, n=43 donors) — B-cell state remodeling and clonal dynamics
3. **Lung spatial transcriptomics + BALF scRNA-seq** (continuity analysis from GSE223793/PRJCA042363) — FOXO3-driven B-cell niches in injured lung

## Repository Structure

```
.
├── README.md
├── 01_proteomics/
│   ├── 01_QC_protein_count_per_sample.R
│   ├── 02_PCA_visualization.R
│   ├── 03_differential_protein_analysis.R
│   ├── 04_heatmap_control_DEPs.R
│   ├── 05_heatmap_immunoglobulin_comparison.R
│   ├── 06_WGCNA_module_analysis.R
│   ├── 07_ARDS_protein_expression_boxplot.R
│   └── 08_prognostic_panel_ROC.R
├── 02_scRNAseq_BCR/
│   ├── 01_cell_proportion_calculation.R
│   ├── 02_BCR_clonotype_matrix_extraction.R
│   ├── 03_GO_enrichment_bubble_plot.R
│   └── 04_SCENIC_transcription_factor_analysis.R
└── 03_spatial_transcriptomics/
    ├── 01_gene_set_scoring_spatial.R
    └── 02_spatial_Bcell_niche_analysis.R
```

## Data Availability

| Dataset | Accession | Description |
|---------|-----------|-------------|
| Lung scRNA-seq (prior study) | [GSE223793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223793) | Pediatric ARDS lung multi-omics atlas |
| Spatial transcriptomics | [PRJCA042363](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA042363) | Visium FFPE spatial data |
| Peripheral blood scRNA-seq + BCR | Pending | Influenza-associated ARDS cohort |
| Plasma proteomics | Pending | DIA-MS quantification matrix |

## Citation

> Song L, et al. Age-dependent FOXO3 rewiring of B-cell immunity shapes outcomes in viral acute respiratory distress syndrome. *(Under review, 2025)*

## License

MIT License
