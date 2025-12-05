import pandas
from google.cloud import storage


bucket_name = 'pcst-ad'
client = storage.Client()
bucket = client.bucket(bucket_name)

# list celltypes by finding network.tsv files
blobs = list(bucket.list_blobs(prefix='output/'))
network_blobs = [b for b in blobs if b.name.endswith('/network.tsv')]

stats = []

for blob in network_blobs:
    celltype = blob.name.split('/')[1]
    
    # read directly from GCS
    content = blob.download_as_text()
    df = pandas.read_csv(pandas.io.common.StringIO(content), sep='\t')
    
    fava_scores = df['fava_score']
    
    stats.append({
        'celltype': celltype,
        'edge_count': len(df),
        'fava_min': fava_scores.min(),
        'fava_max': fava_scores.max(),
        'fava_mean': fava_scores.mean(),
        'fava_median': fava_scores.median(),
        'fava_std': fava_scores.std(),
        'fava_q25': fava_scores.quantile(0.25),
        'fava_q75': fava_scores.quantile(0.75)
    })

summary = pandas.DataFrame(stats)
print(summary)

summary.to_csv('network_summary.tsv', sep='\t', index=False)

