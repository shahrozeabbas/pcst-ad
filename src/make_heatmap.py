import os
import pickle
import pandas
import seaborn
import matplotlib.pyplot as plt


def jaccard(a, b):
    union = len(a | b)
    intersection = len(a & b)
    return intersection / union if union > 0 else 0


modules = {}
models = snakemake.input.models
celltypes = snakemake.params.celltypes

for celltype, model_path in zip(celltypes, models):
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    modules[celltype] = set(model.module.vs['name'])

# build similarity matrix
matrix = pandas.DataFrame(index=celltypes, columns=celltypes, dtype=float)

for ct1 in celltypes:
    for ct2 in celltypes:
        matrix.loc[ct1, ct2] = jaccard(modules[ct1], modules[ct2])


# plot heatmap
plt.figure(figsize=(8, 6))

seaborn.heatmap(
    matrix.astype(float),
    annot=True,
    fmt='.2f',
    cmap='Blues',
    vmin=0,
    vmax=1,
    square=True
)

plt.title('Module Gene Jaccard Similarity')
plt.tight_layout()

plt.savefig(snakemake.output.heatmap, dpi=150)
