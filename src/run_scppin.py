import pandas
import scppin
import pickle


method = snakemake.params.method
network = pandas.read_csv(snakemake.input.network, sep='\t')


if method == 'mean':
    network['weight'] = (network['fava_score'] + network['string_score']) / 2
elif method == 'multiply':
    network['weight'] = (network['fava_score'] * network['string_score'])
elif method == 'geometric_mean':
    network['weight'] = (network['fava_score'] * network['string_score']) ** 0.5
elif method == 'fava':
    network['weight'] = network['fava_score']
elif method == 'string':
    network['weight'] = network['string_score']
else:
    raise ValueError(f'Invalid method: {method}')


pvals_df = pandas.read_csv(snakemake.input.pvalues, sep='\t', index_col=0)
pvalues = dict(zip(pvals_df['GENE'], pvals_df['P']))


model = scppin.scPPIN()

model.load_network(network, weight_column='weight')

model.set_node_weights(pvalues)

model.detect_module(
    fdr=0.01, edge_weight_attr='weight', 
    normalization=None, use_max_prize_root=True
)

with open(snakemake.output.model, 'wb') as f:
    pickle.dump(model, f)

print(f'Nodes: {model.module.vcount()}')
print(f'Edges: {model.module.ecount()}')