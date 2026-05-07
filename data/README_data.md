# Data Sources and Input File Descriptions

This file describes all input data required to run the analysis scripts in this repository. **Raw data files are not included in this repository** due to patient privacy considerations and file size constraints. Processed and anonymized datasets will be deposited in public repositories prior to publication.

---

## 1. Plasma Proteomics Input

| File | Format | Description |
|------|--------|-------------|
| `ALL_sample.final_matrix.csv` | CSV | FOT-normalized protein abundance matrix. Rows = protein/gene symbols (HGNC); Columns = sample IDs (S1–S108, BF1–BF5). Values represent FOT × 10⁵. |
| `ALL_sample.raw_matrix.csv` | CSV | Raw intensity matrix before normalization (used for QC). |

**Sample ID to group mapping:**

| Sample IDs | Group code | Clinical group |
|------------|------------|----------------|
| S1–S14 | `AC` | Healthy adult controls |
| S15–S27 | `CC` | Pediatric healthy controls |
| S28–S47, S67–S76, S87–S108, BF1–BF5 | `AA` | Adult ARDS (AARDS) |
| S48–S57, S77–S86 | `CA1` | Pediatric ARDS — acute phase (T1) |
| S58–S66 | `CA2` | Pediatric ARDS — outcome phase (T2) |

---

## 2. Single-Cell RNA-seq and BCR Input

| Object | Format | Description |
|--------|--------|-------------|
| `CountsSeuratB_umap_J` | Seurat RDS | B-cell-focused Seurat object after Harmony integration. Contains metadata columns: `group`, `group1`, `celltype2`, `sample`, `orig.ident`. |
| `meta_data` | data.frame (CSV or in-memory) | Full cell metadata for all lineages. Columns: `orig.ident`, `group`, `celltype_final`, `celltype_fine`. |
| `scBCR_RNA_PB` | Seurat/scRepertoire object | BCR-annotated Seurat object. Contains: `cloneType`, `Isotype`, `sample`, `group1`. |

**Single-cell group label conventions:**

| Metadata label | Clinical group |
|----------------|----------------|
| `childcontrol` / `CC` | Pediatric healthy controls |
| `adultcontrol` / `AC` | Healthy adult controls |
| `childARDS_T1` / `CA` | Pediatric ARDS — acute phase |
| `adultARDS_T1` / `AA` | Adult ARDS — acute phase |

---

## 3. Spatial Transcriptomics Input

| Object | Format | Description |
|--------|--------|-------------|
| `Spatial_integrated` | Seurat RDS | Integrated Visium FFPE Seurat object from prior publication. Contains assays: `RNA` (or `SCT`), `Spatial`, `celltypeprops` (NNLS/RCTD deconvolution weights). |
| `Spatial_integrated_Bcell_region_analysis` | Seurat RDS | Subset object used specifically for B-cell niche analysis. |

**Data accession:**

| Dataset | Accession | Access |
|---------|-----------|--------|
| Lung tissue scRNA-seq | [GSE223793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223793) | Public |
| Spatial transcriptomics | [PRJCA042363](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA042363) | Controlled access |

---

## 4. Functional Annotation Input (for GO enrichment scripts)

| File | Source | Description |
|------|--------|-------------|
| `GO_bubble_4.xlsx` | Metascape / clusterProfiler output | Multi-sheet Excel file; each sheet = one cluster's GO enrichment results. Columns: `Description`, `gene_ratio`, `LogP`, `Symbols`. |
| `cisTarget_databases/` | [RcisTarget](https://resources.aertslab.org/cistarget/) | Required for SCENIC motif enrichment. Download separately. |

---

## 5. Deposited Data (upon publication)

Peripheral blood scRNA-seq, BCR sequencing, and plasma proteomics data generated in this study will be deposited in the Gene Expression Omnibus (GEO) and/or the National Genomics Data Center (NGDC) prior to publication. Accession numbers will be updated here.
