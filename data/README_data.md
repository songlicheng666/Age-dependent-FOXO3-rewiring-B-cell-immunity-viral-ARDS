# Data Sources

Input data files are not included in this repository due to patient privacy and file size constraints.

## Proteomics Input
- `ALL_sample.final_matrix.csv`: FOT-normalized protein abundance matrix (rows=proteins, cols=sample IDs S1-S108)
- Groups: AA=Adult ARDS, AC=Adult controls, CA1=Pediatric ARDS T1, CA2=Pediatric ARDS T2, CC=Pediatric controls

## Single-cell Input
- Seurat RDS object: `CountsSeuratB_umap_J` (B-cell focused, Harmony-integrated)
- scRepertoire object: `scBCR_RNA_PB` (BCR clonotype annotations)

## Spatial Transcriptomics Input
- Seurat Visium object: `Spatial_integrated` (from prior publication)
- Accession: GSE223793 (lung scRNA-seq); PRJCA042363 (spatial, controlled access)

## Deposited Data (upon publication)
Peripheral blood scRNA-seq, BCR, and proteomics data will be deposited in GEO/NGDC prior to publication.
