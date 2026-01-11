# pcst-ad

Snakemake pipeline for identifying Alzheimer's disease-relevant protein modules by integrating three data modalities:

- **Edges**: Protein-protein interactions from STRING database
- **Edge weights**: Differential co-expression scores from FAVA (Control vs AD)
- **Node weights**: Gene-level p-values from MAGMA (AD GWAS)

The Prize-Collecting Steiner Tree (PCST) algorithm finds connected subnetworks that balance including high-prize nodes (genetically associated genes) with traversing high-weight edges (differentially co-expressed protein interactions).

## Pipeline Overview

```mermaid
graph LR
    A[counts] --> B[fava]
    B --> C[network]
    C --> D[scppin]
    D --> E[gsea]
    D --> F[heatmap]
```

### Workflow Steps

- **counts**: Extract cell-type specific expression data filtered to MAGMA GWAS genes
- **fava**: Compute gene-gene co-expression scores for Control and AD conditions (GPU-accelerated)
- **network**: Compute differential correlation using Fisher z-transformation, filter edges to known STRING protein-protein interactions
- **scppin**: PCST module detection—edges weighted by differential co-expression, nodes weighted by GWAS significance
- **gsea**: Gene set enrichment analysis on detected modules
- **heatmap**: Cross-cell-type module comparison using Jaccard similarity

## Input Data

Required input files in `input/` directory:

- `magma_gene_symbol_results.tsv` - GWAS gene-level association results from MAGMA
- `GSE174367_snRNA-seq_cell_meta.csv` - Single-cell RNA-seq metadata with cell type and diagnosis annotations
- `GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5` - Expression count matrix in 10x Genomics format

## Configuration

Create a `config.yaml` file:

```yaml
BUCKET: gs://your-bucket

CONDITIONS:
  - Control
  - AD

CELLTYPES:
  - ODC
  - EX
  - INH
  - ASC
  - MG
  - OPC
  - PER.END

# GPU settings for fava rule
GPU_TYPE: nvidia-tesla-t4
GPU_MACHINE_TYPE: n1-highmem-2
GPU_BOOT_IMAGE: projects/rocky-linux-accelerator-cloud/global/images/family/rocky-linux-8-optimized-gcp-nvidia-latest

# GSEA gene sets
GENE_SETS:
  - Reactome_2022
  - DisGeNET
```

## Output Files

### Per Cell Type (`output/{celltype}/`)
- `{condition}_anndata.h5ad` - Filtered AnnData object with highly variable genes
- `{condition}_fava_pairs.tsv` - Gene pair co-expression scores
- `network.tsv` - Differential correlation network filtered by STRING PPI
- `scppin_object.pkl` - PCST module detection results
- `module_gsea_enrichment.tsv` - Pathway enrichment results

### Figures (`figures/`)
- `{celltype}/{condition}_fava_histogram.png` - Distribution of FAVA scores
- `module_jaccard_heatmap.png` - Cross-cell-type module similarity heatmap

## Usage

### Local execution
```bash
snakemake --use-conda
```

### Google Cloud Platform execution
```bash
snakemake --profile profiles/gcp
```
