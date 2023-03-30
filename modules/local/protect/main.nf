process PROTECT {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/protect:2.3--0'

    input:
    tuple val(meta), path(chord_prediction), path(purple_dir), path(linx_dir), path(lilac_dir), path(virusinterpreter)
    val genome_ver
    path serve_resources
    path driver_gene_panel
    path disease_ontology

    output:
    tuple val(meta), path('*.protect.tsv'), emit: tsv
    path 'versions.yml'                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            -tumor_sample_id ${meta.tumor_id} \\
            -reference_sample_id ${meta.normal_id} \\
            -primary_tumor_doids '' \\
            -output_dir ./ \\
            \\
            -serve_actionability_dir ${serve_resources} \\
            -ref_genome_version ${genome_ver} \\
            -driver_gene_tsv ${driver_gene_panel} \\
            -doid_json ${disease_ontology} \\
            \\
            -chord_prediction_txt ${chord_prediction} \\
            \\
            -purple_qc_file ${purple_dir}/${meta.tumor_id}.purple.qc \\
            -purple_purity_tsv ${purple_dir}/${meta.tumor_id}.purple.purity.tsv \\
            -purple_gene_copy_number_tsv ${purple_dir}/${meta.tumor_id}.purple.cnv.gene.tsv \\
            -purple_somatic_driver_catalog_tsv ${purple_dir}/${meta.tumor_id}.driver.catalog.somatic.tsv \\
            -purple_somatic_variant_vcf ${purple_dir}/${meta.tumor_id}.purple.somatic.vcf.gz \\
            -purple_germline_driver_catalog_tsv ${purple_dir}/${meta.tumor_id}.driver.catalog.germline.tsv \\
            -purple_germline_variant_vcf ${purple_dir}/${meta.tumor_id}.purple.germline.vcf.gz \\
            \\
            -linx_breakend_tsv ${linx_dir}/${meta.tumor_id}.linx.breakend.tsv \\
            -linx_driver_catalog_tsv ${linx_dir}/${meta.tumor_id}.linx.driver.catalog.tsv \\
            -linx_fusion_tsv ${linx_dir}/${meta.tumor_id}.linx.fusion.tsv \\
            \\
            -lilac_qc_csv ${lilac_dir}/${meta.tumor_id}.lilac.qc.csv \\
            -lilac_result_csv ${lilac_dir}/${meta.tumor_id}.lilac.csv \\
            \\
            -annotated_virus_tsv ${virusinterpreter}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protect: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*PROTECT v//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_id}.protect.tsv"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
