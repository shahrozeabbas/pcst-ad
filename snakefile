
envvars:
    'CONDA_PLUGINS_AUTO_ACCEPT_TOS'

configfile: 'config.yaml'

CELLTYPES = config['CELLTYPES']
CONDITIONS = config['CONDITIONS']


rule all:
    input:
        expand(
            'output/{celltype}/module_gsea_enrichment.tsv', 
            celltype=CELLTYPES
        ),
        'figures/module_jaccard_heatmap.png',
        expand('figures/{celltype}/module_network.png', celltype=CELLTYPES)
        
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
    resources:
        gpu=1,
        nvidia_gpu=config['GPU_TYPE'],
        googlebatch_machine_type=config['GPU_MACHINE_TYPE'],
        googlebatch_boot_disk_image=config['GPU_BOOT_IMAGE']
    conda:
        'containers/fava/env.yaml'
    script:
        'src/run_fava.py'

rule network:
    input:
        control_edges='output/{celltype}/Control_fava_pairs.tsv',
        disease_edges='output/{celltype}/AD_fava_pairs.tsv'
    output:
        network='output/{celltype}/network.tsv'
    params:
        cutoff=config['STRING_CUTOFF']
    resources:
        googlebatch_memory=32000,
        googlebatch_machine_type='e2-highmem-4'
    conda:
        'containers/network/env.yaml'
    script:
        'src/make_network.py'

rule scppin:
    input:
        network='output/{celltype}/network.tsv',
        pvalues='input/magma_gene_symbol_results.tsv'
    output:
        model='output/{celltype}/scppin_object.pkl'
    params:
        method='fava'
    resources:
        googlebatch_memory=32000,
        googlebatch_machine_type='e2-highmem-4'
    conda:
        'containers/scppin/env.yaml'
    script:
        'src/run_scppin.py'

rule plot:
    input:
        model='output/{celltype}/scppin_object.pkl',
        pvalues='input/magma_gene_symbol_results.tsv'
    output:
        plot='figures/{celltype}/module_network.png'
    resources:
        googlebatch_memory=32000,
        googlebatch_machine_type='e2-highmem-4'
    conda:
        'containers/plot/env.yaml'
    script:
        'src/make_plot.py'

rule gsea:
    input:
        model='output/{celltype}/scppin_object.pkl'
    output:
        enrichment='output/{celltype}/module_gsea_enrichment.tsv'
    params:
        gene_sets=config['GENE_SETS']
    conda:
        'containers/gsea/env.yaml'
    script:
        'src/run_gsea.py'

rule heatmap:
    input:
        models=expand('output/{celltype}/scppin_object.pkl', celltype=CELLTYPES)
    output:
        heatmap='figures/module_jaccard_heatmap.png'
    params:
        celltypes=CELLTYPES
    resources:
        googlebatch_memory=32000,
        googlebatch_machine_type='e2-highmem-4'
    conda:
        'containers/heatmap/env.yaml'
    script:
        'src/make_heatmap.py'
