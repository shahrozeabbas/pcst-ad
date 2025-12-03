
envvars:
    'CONDA_PLUGINS_AUTO_ACCEPT_TOS'

configfile: 'config.yaml'

CELLTYPES = config['CELLTYPES']
CONDITIONS = config['CONDITIONS']


rule all:
    input:
        expand(
            'output/{celltype}/{condition}_fava_pairs.tsv', 
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
    script:
        'src/get_counts.py'

rule fava:
    input:
        anndata='output/{celltype}/{condition}_anndata.h5ad'
    output:
        pairs='output/{celltype}/{condition}_fava_pairs.tsv',
        histogram='figures/{celltype}/{condition}_fava_histogram.png'
    params:
        cc_cutoff=0.0
    resources:
        nvidia_gpu=config['GPU_TYPE'],
        googlebatch_machine_type=config['GPU_MACHINE_TYPE'],
        googlebatch_boot_disk_image=config['GPU_BOOT_IMAGE']
    conda:
        'containers/fava/env.yaml'
    script:
        'src/run_fava.py'

rule merge:
    input:
        first_edges='output/{celltype}/Control_fava_pairs.tsv',
        second_edges='output/{celltype}/Disease_fava_pairs.tsv'
    output:
        edges='output/{celltype}/merged_fava.tsv'
    conda:
        'containers/merge/polars.yaml'
    script:
        'src/merge_fava.py'

rule stringdb:
    input:
        pairs='output/{celltype}/merged_fava.tsv'
    output:
        network='output/{celltype}/network.tsv'
    conda:
        'containers/stringdb/env.yaml'
    script:
        'src/make_network.py'
