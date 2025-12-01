import polars as pl


def prep(df, a, b):
    return df.with_columns([
        pl.min_horizontal(a, b).alias('u'),
        pl.max_horizontal(a, b).alias('v')
    ])


# read TSV edge tables
df1 = pl.read_csv(snakemake.input.first_edges, sep='\t')
df2 = pl.read_csv(snakemake.input.second_edges, sep='\t')

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
        (pl.col('score1') - pl.col('score2')).abs().alias('merge_score')
    ])
)

edges.write_csv(snakemake.output.edges, sep='\t')
