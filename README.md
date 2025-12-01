# pcst-ad

Snakemake pipeline for PCST (Prize-Collecting Steiner Tree) analysis of Alzheimer's disease single-cell RNA-seq data. Identifies protein-protein interaction networks using FAVA and STRING database.

## Configuration

Create a `config.yaml` file:

```yaml
BUCKET: gs://your-bucket

CONDITIONS:
  - Control
  - Disease

CELLTYPES:
  - Neuron
  - Astrocyte
  - Microglia
```

## Usage

```bash
snakemake --use-conda
```

