// TODO(SW): accept multiple BAMs
// TODO(SW): appropriately select preferential input source (can be either DNA or RNA BAM); manually set this in channel to processing?
// TODO(SW): gather labels, BAMs, and BAIs into single entry for SAGE append
// TODO(SW): adjust meta_append to have appropriate IDs
// TODO(SW): update markdups for wgts and targeted workflows




//
// SAGE append adds additional sample data to an existing SAGE VCF
//

import Constants
import Utils

include { SAGE_APPEND as SOMATIC } from '../../../modules/local/sage/append/main'
include { SAGE_APPEND as GERMLINE } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND {
    take:
    // Sample data
    ch_inputs     // channel: [mandatory] [ meta ]
    ch_purple_dir // channel: [mandatory] [ meta, purple_dir ]
    ch_bam        // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_version   // channel: [mandatory] genome version
    genome_fai       // channel: [mandatory] /path/to/genome_fai
    genome_dict      // channel: [mandatory] /path/to/genome_dict

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_bam,
        ch_purple_dir,
    )
        .map { meta, bam, bai, purple_dir ->
            return [
                meta,



                // TODO(SW): solve this

                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Utils.selectCurrentOrExisting(bai, meta, Constants.INPUT.BAI_RNA_TUMOR),



                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, bam, bai, purple_dir ->
            runnable: bam && purple_dir
            skip: true
                return meta
        }

    //
    // MODULE: SAGE append germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, bam, bai, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, bam, bai, purple_dir ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta)

            def has_normal_dna = Utils.hasNormalDna(meta)




            // TODO(SW): solve this and below

            def has_tumor_rna = Utils.hasTumorRna(meta)




            def has_smlv_germline = file(purple_dir).resolve("${tumor_dna_id}.purple.germline.vcf.gz")
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_VCF_NORMAL)

            runnable: has_normal_dna && has_tumor_rna && has_smlv_germline && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, bam, bai ]
    ch_sage_append_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, bam, bai, purple_dir ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta)

            def meta_append = [
                key: meta.group_id,
                id: meta.group_id,




                // TODO(SW): solve this

                tumor_rna_id: Utils.getTumorRnaSampleName(meta),




                dna_id: Utils.getNormalDnaSampleName(meta),
            ]

            def purple_smlv_vcf = file(purple_dir).resolve("${tumor_dna_id}.purple.germline.vcf.gz")

            return [meta_append, purple_smlv_vcf, bam, bai]
        }

    // Run process
    GERMLINE(
        ch_sage_append_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
    )

    ch_versions = ch_versions.mix(GERMLINE.out.versions)

    //
    // MODULE: SAGE append somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, bam, bai, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, bam, bai, purple_dir ->
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta)

            def has_tumor_dna = Utils.hasTumorDna(meta)





            // TODO(SW): solve this

            def has_tumor_rna = Utils.hasTumorRna(meta)





            def has_smlv_somatic = file(purple_dir).resolve("${tumor_dna_id}.purple.somatic.vcf.gz")
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_VCF_TUMOR)

            runnable: has_tumor_dna && has_tumor_rna && has_smlv_somatic && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, bam, bai]
    ch_sage_append_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, bam, bai, purple_dir ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta)

            def meta_append = [
                key: meta.group_id,
                id: meta.group_id,




                // TODO(SW): solve this

                tumor_rna_id: Utils.getTumorRnaSampleName(meta),





                dna_id: Utils.getTumorDnaSampleName(meta),
            ]

            def purple_smlv_vcf = file(purple_dir).resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            return [meta_append, purple_smlv_vcf, bam, bai]
        }

    // Run process
    SOMATIC(
        ch_sage_append_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
    )

    ch_versions = ch_versions.mix(SOMATIC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_vcf ]
    ch_somatic_vcf = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    ch_germline_vcf = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic_vcf  = ch_somatic_vcf  // channel: [ meta, sage_append_vcf ]
    germline_vcf = ch_germline_vcf // channel: [ meta, sage_append_vcf ]

    versions     = ch_versions     // channel: [ versions.yml ]
}
