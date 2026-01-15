import pickle
import math
import random
import tempfile
import pandas
import scppin
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from PIL import Image
from adjustText import adjust_text

CELLTYPE_NAMES = {
    'ODC': 'Oligodendrocytes',
    'EX': 'Excitatory Neurons',
    'INH': 'Inhibitory Neurons',
    'ASC': 'Astrocytes',
    'MG': 'Microglia',
    'OPC': 'Oligodendrocyte Precursor Cells',
    'PER.END': 'Pericytes/Endothelial'
}

celltype = snakemake.wildcards.celltype
celltype_name = CELLTYPE_NAMES.get(celltype, celltype).upper()

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
colors_hex = ['#{:02x}{:02x}{:02x}'.format(
    int(r * 255), int(grn * 255), int(b * 255)
) for r, grn, b, _ in colors]

random.seed(0)
layout = g.layout('lgl')
layout.scale(1.5)
coords = layout.coords

visual_style = {
    'layout': layout,
    'vertex_size': 30,
    'vertex_color': colors_hex,
    'vertex_frame_color': '#2C5282',
    'vertex_frame_width': 1,
    'edge_color': '#A0AEC0',
    'edge_width': 0.8,
    'edge_curved': 0.1,
    'margin': 80
}

with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
    tmp_path = tmp.name

ig.plot(g, tmp_path, bbox=(1200, 700), **visual_style)
network_img = Image.open(tmp_path)

coords_x = [c[0] for c in coords]
coords_y = [c[1] for c in coords]
x_min, x_max = min(coords_x), max(coords_x)
y_min, y_max = min(coords_y), max(coords_y)

x_range = x_max - x_min
y_range = y_max - y_min
x_padding = x_range * 0.1
y_padding = y_range * 0.1
center_x = sum(coords_x) / len(coords_x)
center_y = sum(coords_y) / len(coords_y)
offset_scale = 0.02 * max(x_range, y_range)

fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 1, height_ratios=[20, 1], hspace=0.3)
ax = fig.add_subplot(gs[0])
cbar_ax = fig.add_subplot(gs[1])

ax.imshow(network_img, extent=[x_min - x_padding, x_max + x_padding,
                                y_max + y_padding, y_min - y_padding],
          aspect='auto', zorder=0)

texts = []
for (x, y), name in zip(coords, gene_names):
    dx = x - center_x
    dy = y - center_y
    norm = math.hypot(dx, dy)
    if norm > 0:
        x += (dx / norm) * offset_scale
        y += (dy / norm) * offset_scale
    texts.append(
        ax.text(x, y, name, fontsize=10, color='#1A202C', weight='bold', zorder=3)
    )
adjust_text(
    texts,
    ax=ax,
    x=coords_x,
    y=coords_y,
    expand_text=(1.2, 1.35),
    expand_points=(1.5, 1.5),
    force_text=(0.18, 0.18),
    force_points=(0.35, 0.35),
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
)

ax.set_xlim(x_min - x_padding, x_max + x_padding)
ax.set_ylim(y_min - y_padding, y_max + y_padding)
ax.axis('off')
ax.set_title(celltype_name, fontsize=24, weight='bold', pad=20, loc='left')

cmap = cm.Blues
norm = mcolors.Normalize(vmin=0, vmax=max_neg_log)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('-log10(p-value)', fontsize=10)

plt.savefig(snakemake.output.plot, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
