import polars
import pandas
import stringdb


def canonicalize(df, a, b):
    return df.with_columns([
        polars.min_horizontal(a, b).alias('u'),
        polars.max_horizontal(a, b).alias('v')
    ])


# read both condition edge tables
df1 = polars.read_csv(snakemake.input.control_edges, separator='\t')
df2 = polars.read_csv(snakemake.input.disease_edges, separator='\t')

# canonicalize undirected edges
df1 = canonicalize(df1, 'Protein_1', 'Protein_2')
df2 = canonicalize(df2, 'Protein_1', 'Protein_2')

# intersect fava edges across conditions, compute |control - disease|
fava_edges = (
    df1.join(
        df2.select(['u', 'v', 'Score']).rename({'Score': 'disease_score'}),
        on=['u', 'v'],
        how='inner'
    )
    .rename({'Score': 'control_score'})
    .with_columns([
        (polars.col('control_score') - polars.col('disease_score')).abs().alias('fava_score')
    ])
)

# get unique proteins for stringdb query
proteins = list(set(fava_edges['u'].to_list() + fava_edges['v'].to_list()))

# query stringdb in batches
string_ids_df = stringdb.get_string_ids(proteins, species=9606)
string_ids = string_ids_df['stringId'].dropna().drop_duplicates().tolist()

edges = []
batch_size = 1000

for i in range(0, len(string_ids), batch_size):
    batch = string_ids[i:i + batch_size]
    print(f'Processing batch {i // batch_size + 1}: proteins {i}-{min(i + batch_size, len(string_ids))}')

    string_edges = stringdb.get_interaction_partners(
        identifiers=batch, species=9606,
        required_score=400, limit=None
    )

    print(f'  Got {len(string_edges)} interactions in this batch')
    edges.append(string_edges)

if edges:
    string_df = polars.from_pandas(pandas.concat(edges, ignore_index=True))
    string_df = canonicalize(string_df, 'preferredName_A', 'preferredName_B')
    string_df = string_df.select(['u', 'v', 'score']).unique(subset=['u', 'v'])

    # final intersection: fava edges that exist in stringdb
    network = fava_edges.join(
        string_df.select(['u', 'v', 'score']).rename({'score': 'string_score'}),
        on=['u', 'v'],
        how='inner'
    )

    network.select(['u', 'v', 'fava_score', 'string_score']).write_csv(
        snakemake.output.network, separator='\t'
    )
else:
    print('No interactions retrieved')
