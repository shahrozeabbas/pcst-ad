
configfile: 'config.yaml'


BUCKET = config['BUCKET']

CELLTYPES = config['CELLTYPES']
CONDITIONS = config['CONDITIONS']


rule all:
    input:
        expand(
            f'{BUCKET}/output/{{celltype}}/network.tsv',
            celltype=CELLTYPES
        ),
        expand(
            f'{BUCKET}/figures/{{celltype}}/{{condition}}_fava_histogram.png',
            celltype=CELLTYPES,
            condition=CONDITIONS
        )

rule counts:
    input:
        gwas=f'{BUCKET}/input/magma_gene_symbol_results.tsv',
        metadata=f'{BUCKET}/input/GSE174367_snRNA-seq_cell_meta.csv',
        counts=f'{BUCKET}/input/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5'
    output:
        anndata=f'{BUCKET}/output/{{celltype}}/{{condition}}_anndata.h5ad'
    conda:
        'envs/scanpy.yaml'
    script:
        'src/get_counts.py'

rule fava:
    input:
        anndata=f'{BUCKET}/output/{{celltype}}/{{condition}}_anndata.h5ad'
    output:
        pairs=f'{BUCKET}/output/{{celltype}}/{{condition}}_fava_pairs.tsv',
        histogram=f'{BUCKET}/figures/{{celltype}}/{{condition}}_fava_histogram.png'
    params:
        cc_cutoff=0.3
    conda:
        'envs/fava.yaml'
    script:
        'src/run_fava.py'

rule merge_fava:
    input:
        first_edges=f'{BUCKET}/output/{{celltype}}/Control_fava_pairs.tsv',
        second_edges=f'{BUCKET}/output/{{celltype}}/Disease_fava_pairs.tsv'
    output:
        edges=f'{BUCKET}/output/{{celltype}}/merged_fava.tsv'
    conda:
        'envs/polars.yaml'
    script:
        'src/merge_fava.py'

rule stringdb:
    input:
        pairs=f'{BUCKET}/output/{{celltype}}/merged_fava.tsv'
    output:
        network=f'{BUCKET}/output/{{celltype}}/network.tsv'
    conda:
        'envs/stringdb.yaml'
    script:
        'src/make_network.py'
