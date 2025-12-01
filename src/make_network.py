import numpy
import pandas
import stringdb


pairs = pandas.read_csv(snakemake.input.pairs, sep='\t')

proteins = list(
    set(
        pairs['Protein_1'].tolist()
        + pairs['Protein_2'].tolist()
    )
)


string_ids_df = stringdb.get_string_ids(proteins, species=9606)

string_ids = (
    string_ids_df['stringId']
    .dropna()
    .drop_duplicates()
    .tolist()
)


edges = []
batch_size = 1000

for i in range(0, len(string_ids), batch_size):

    batch = string_ids[i:i + batch_size]
    print(f'Processing batch {i // batch_size + 1}: proteins {i}-{min(i + batch_size, len(string_ids))}')

    string_edges = stringdb.get_interaction_partners(
        identifiers=batch, species=9606,
        required_score=400, limit=None
    )

    print('  Got', len(string_edges), 'interactions in this batch')
    edges.append(string_edges)


if edges:

    string_df = pandas.concat(edges, ignore_index=True)
    string_df = string_df[['preferredName_A', 'preferredName_B', 'score']].copy()
    string_df.columns = ['protein_1', 'protein_2', 'confidence_score']

    string_df = string_df.drop_duplicates(subset=['protein_1', 'protein_2'], keep='first')
    string_df = string_df.sort_values('confidence_score', ascending=False).reset_index(drop=True)

    df1 = pairs.rename(columns={'Protein_1':'a','Protein_2':'b','Score':'fava_score'})
    df2 = string_df.rename(columns={'protein_1':'a','protein_2':'b','confidence_score':'string_score'})

    df1[['u','v']] = numpy.sort(df1[['a','b']], axis=1)
    df2[['u','v']] = numpy.sort(df2[['a','b']], axis=1)

    network = df1.merge(
        df2[['u','v','string_score']],
        on=['u','v'],
        how='inner'
    )

    network = network[['u','v','fava_score','string_score']]
    network.to_csv(snakemake.output.network, sep='\t', index=False)

else:
    print('No interactions retrieved')


