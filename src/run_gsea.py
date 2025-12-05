import pandas
import gseapy
import pickle


with open(snakemake.input.model, 'rb') as f:
    model = pickle.load(f)

module_genes = model.module.vs['name']
background_genes = model.network.vs['name']

enr = gseapy.enrichr(
    gene_list=list(module_genes),
    gene_sets=snakemake.params.gene_sets,
    background=list(background_genes),
    organism='human',
    outdir=None
)

enr.results.to_csv(snakemake.output.enrichment, sep='\t', index=False)
