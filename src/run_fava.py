import scanpy
import pandas
from favapy import fava
import matplotlib.pyplot as plt


adata = scanpy.read_h5ad(snakemake.input.anndata)


pairs = fava.cook(
    data=adata, log2_normalization=False, 
    epochs=20, batch_size=32, random_seed=0,
    CC_cutoff=snakemake.params.cc_cutoff, layer='counts'
)

pairs['Score'].plot(
    kind='hist', bins=50, edgecolor='black', 
    alpha=0.8, title='Distribution of Fava Score'
)

plt.xlabel('Score')
plt.ylabel('Frequency')


pairs.to_csv(snakemake.output.pairs, sep='\t', index=False)
plt.savefig(snakemake.output.histogram, dpi=300, bbox_inches='tight')
