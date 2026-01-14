import pickle
import pandas
import math
import io
import tempfile
import scppin
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from PIL import Image, ImageDraw, ImageFont


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
# Convert RGBA (0-1) to hex strings for igraph
colors = ['#{:02x}{:02x}{:02x}'.format(
    int(r * 255), int(grn * 255), int(b * 255)
) for r, grn, b, _ in colors]

layout = g.layout('kamada_kawai')
layout.scale(1.5)  # Spread nodes apart to reduce label overlap

# Compute label angles for radial placement (labels point outward from center)
layout_coords = layout.coords
center_x = sum(c[0] for c in layout_coords) / len(layout_coords)
center_y = sum(c[1] for c in layout_coords) / len(layout_coords)
label_angles = [math.atan2(y - center_y, x - center_x) for x, y in layout_coords]

visual_style = {
    'layout': layout,
    'vertex_size': 25,
    'vertex_color': colors,
    'vertex_frame_color': '#2C5282',
    'vertex_frame_width': 1,
    'vertex_label': gene_names,
    'vertex_label_size': 15,
    'vertex_label_color': '#1A202C',
    'vertex_label_dist': 1.5,
    'vertex_label_angle': label_angles,
    'edge_color': '#A0AEC0',
    'edge_width': 0.8,
    'edge_curved': 0.1,
    'margin': 80
}

# Render network with Cairo to temp file
with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
    tmp_path = tmp.name
ig.plot(g, tmp_path, bbox=(1200, 700), **visual_style)
network_img = Image.open(tmp_path)

# Create colorbar with matplotlib
fig, ax = plt.subplots(figsize=(10, 0.5))
cmap = cm.Blues
norm = mcolors.Normalize(vmin=0, vmax=max_neg_log)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, cax=ax, orientation='horizontal')
cbar.set_label('-log10(p-value)', fontsize=10)
plt.tight_layout()

# Save colorbar to buffer
cbar_buffer = io.BytesIO()
plt.savefig(cbar_buffer, format='png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
cbar_buffer.seek(0)
cbar_img = Image.open(cbar_buffer)

# Resize colorbar to match network width
cbar_img = cbar_img.resize((network_img.width, int(cbar_img.height * network_img.width / cbar_img.width)))

# Composite images
final_height = network_img.height + cbar_img.height + 20
final_img = Image.new('RGB', (network_img.width, final_height), 'white')
final_img.paste(network_img, (0, 0))
final_img.paste(cbar_img, (0, network_img.height + 20))

# Add cell type name in top left
draw = ImageDraw.Draw(final_img)
try:
    font = ImageFont.truetype('DejaVuSans-Bold.ttf', 24)
except OSError:
    font = ImageFont.load_default()
draw.text((20, 15), celltype_name, fill='black', font=font)

final_img.save(snakemake.output.plot, dpi=(150, 150))
