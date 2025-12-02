
configfile: 'config.yaml'


CELLTYPES = config['CELLTYPES']
CONDITIONS = config['CONDITIONS']


rule all:
    input:
        expand(
            'output/{celltype}/network.tsv',
            celltype=CELLTYPES
        ),
        expand(
            'figures/{celltype}/{condition}_fava_histogram.png',
            celltype=CELLTYPES,
            condition=CONDITIONS
        )

rule counts:
    input:
        gwas='input/magma_gene_symbol_results.tsv',
        metadata='input/GSE174367_snRNA-seq_cell_meta.csv',
        counts='input/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5'
    output:
        anndata='output/{celltype}/{condition}_anndata.h5ad'
    conda:
        'containers/counts/env.yaml'
    resources:
        mem_mb=16000, cpus=2, disk_mb=8000
    script:
        'src/get_counts.py'

rule fava:
    input:
        anndata='output/{celltype}/{condition}_anndata.h5ad'
    output:
        pairs='output/{celltype}/{condition}_fava_pairs.tsv',
        histogram='figures/{celltype}/{condition}_fava_histogram.png'
    params:
        cc_cutoff=0.3
    conda:
        'envs/fava.yaml'
    script:
        'src/run_fava.py'

rule merge:
    input:
        first_edges='output/{celltype}/Control_fava_pairs.tsv',
        second_edges='output/{celltype}/Disease_fava_pairs.tsv'
    output:
        edges='output/{celltype}/merged_fava.tsv'
    conda:
        'envs/polars.yaml'
    script:
        'src/merge_fava.py'

rule stringdb:
    input:
        pairs='output/{celltype}/merged_fava.tsv'
    output:
        network='output/{celltype}/network.tsv'
    conda:
        'envs/stringdb.yaml'
    script:
        'src/make_network.py'
