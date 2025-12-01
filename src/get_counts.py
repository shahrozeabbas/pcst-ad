import scanpy
import pandas


group = snakemake.wildcards.group
celltype = snakemake.wildcards.celltype
adata = scanpy.read_h5ad(snakemake.input.counts)
metadata = pandas.read_csv(snakemake.input.metadata, sep='\t', index_col=0)


adata.var_names_make_unique()
metadata = metadata.set_index('Barcode')
adata = adata[adata.obs_names.isin(metadata.index)].copy()
adata.obs = adata.obs.join(metadata, how='left') 


gwas = pandas.read_csv(snakemake.input.gwas, sep='\t', index_col=0)

adata = adata[:, adata.var_names.isin(gwas['GENE'])].copy()


mask = (adata.obs['Cell.Type'] == celltype) & (adata.obs['Diagnosis'] == group)

adata = adata[mask].copy()


adata.layers['counts'] = adata.X.copy() # type: ignore

scanpy.pp.highly_variable_genes(
    adata, batch_key='SampleID', layer='counts', 
    flavor='seurat_v3', n_top_genes=10000, subset=True
)


adata.write_h5ad(snakemake.output.anndata)