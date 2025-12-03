import polars


def prep(df, a, b):
    return df.with_columns([
        polars.min_horizontal(a, b).alias('u'),
        polars.max_horizontal(a, b).alias('v')
    ])


# read TSV edge tables
df1 = polars.read_csv(snakemake.input.first_edges, sep='\t')
df2 = polars.read_csv(snakemake.input.second_edges, sep='\t')

# canonicalize undirected endpoints
df1_p = prep(df1, 'protein_1', 'protein_2')
df2_p = prep(df2, 'protein_1', 'protein_2')

# intersect edges on (u, v) and compute |score1 - score2|
edges = (
    df1_p.join(
        df2_p.select(['u', 'v', 'score2']),
        on=['u', 'v'],
        how='inner'
    )
    .select(['u', 'v', 'score1', 'score2'])
    .with_columns([
        (polars.col('score1') - polars.col('score2')).abs().alias('merge_score')
    ])
)

edges.write_csv(snakemake.output.edges, sep='\t')
