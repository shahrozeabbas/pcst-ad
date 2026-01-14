import pickle
import pandas
import math
import scppin
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.cm as cm


with open(snakemake.input.model, 'rb') as f:
    model = pickle.load(f)

g = model.module

if g.vcount() == 0:
    raise ValueError('Module graph is empty (no nodes)')

pvals_df = pandas.read_csv(snakemake.input.pvalues, sep='\t', index_col=0)
pvalues_dict = dict(zip(pvals_df['GENE'], pvals_df['P']))

gene_names = g.vs['name']
node_pvalues = [pvalues_dict.get(gene, 1.0) for gene in gene_names]

neg_log_pvalues = [-math.log10(p) if p > 0 else 0 for p in node_pvalues]

max_neg_log = max(neg_log_pvalues) if neg_log_pvalues else 1.0
if max_neg_log == 0:
    max_neg_log = 1.0
norm_pvalues = [p / max_neg_log if max_neg_log > 0 else 0 for p in neg_log_pvalues]

colors = [cm.Blues(0.3 + 0.7 * p) for p in norm_pvalues]
# Convert RGBA (0-1) to hex strings for igraph
colors = ['#{:02x}{:02x}{:02x}'.format(
    int(r * 255), int(grn * 255), int(b * 255)
) for r, grn, b, _ in colors]

layout = g.layout('fruchterman_reingold', niter=1000)

visual_style = {
    'layout': layout,
    'vertex_size': 25,
    'vertex_color': colors,
    'vertex_frame_color': '#2C5282',
    'vertex_frame_width': 1,
    'vertex_label': gene_names,
    'vertex_label_size': 8,
    'vertex_label_color': '#1A202C',
    'edge_color': '#A0AEC0',
    'edge_width': 0.8,
    'edge_curved': 0.1,
    'margin': 50
}

ig.plot(g, snakemake.output.plot, bbox=(800, 800), **visual_style)
