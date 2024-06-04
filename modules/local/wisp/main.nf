process WISP {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/wisp:1.1.rc1--1'

    input:
    tuple val(meta), path(vcf), path(cobalt_dir), path(purple_dir)
    path genome_fasta
    path genome_fai

    output:
    path 'wisp'
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp /opt/wisp/wisp.jar com.hartwig.hmftools.wisp.purity.PurityEstimator \\
            ${args} \\
            -patient_id ${meta.patient_id} \\
            -tumor_id ${meta.primary_tumor_id} \\
            -ctdna_samples ${meta.sample_id} \\
            \\
            -somatic_vcf ${vcf} \\
            -cobalt_dir ${cobalt_dir} \\
            -purple_dir ${purple_dir} \\
            \\
            -ref_genome ${genome_fasta} \\
            \\
            -gc_ratio_min 0 \\
            -write_types ALL \\
            -log_debug \\
            \\
            -output_dir wisp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisp: \$(java -jar /opt/wisp/wisp.jar -version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p wisp/
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.cn_plot_calcs.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.cn_segments.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.somatic_peak.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.somatic_variants.tsv
    touch wisp/${meta.patient_id}_${meta.sample_id}.wisp.summary.tsv
    touch wisp/${meta.sample_id}.cn_gc_ratio_fit.png
    touch wisp/${meta.sample_id}.somatic_vaf.png

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
