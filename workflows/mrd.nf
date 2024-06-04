import Constants
import Processes
import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse input samplesheet
// NOTE(SW): this is done early and outside of gpars so that we can access synchronously and prior to pipeline execution
inputs = Utils.parseInput(params.input, workflow.stubRun, log)

// Get run config
run_config = WorkflowMain.getRunConfig(params, inputs, log)

// Validate inputs
Utils.validateInput(inputs, run_config, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
]

// TODO(SW): consider whether we should check for null entries here for errors to be more informative
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { AMBER_PROFILING    } from '../subworkflows/local/amber_profiling'
include { COBALT_PROFILING   } from '../subworkflows/local/cobalt_profiling'
include { PREPARE_REFERENCE  } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA } from '../subworkflows/local/read_alignment_dna'
include { READ_PROCESSING    } from '../subworkflows/local/read_processing'
include { SAGE_APPEND        } from '../subworkflows/local/sage_append'
include { WISP_ANALYSIS      } from '../subworkflows/local/wisp_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)

workflow MRD {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = Channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(
        PREPARE_REFERENCE.out.versions,
    )

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_align_dna_tumor_out = Channel.empty()
    if (run_config.stages.alignment) {

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records,
        )

        ch_versions = ch_versions.mix(
            READ_ALIGNMENT_DNA.out.versions,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.dna_tumor)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }

    }

    //
    // SUBWORKFLOW: Run MarkDups for DNA BAMs
    //
    // channel: [ meta, bam, bai ]
    ch_process_dna_tumor_out = Channel.empty()
    if (run_config.stages.markdups) {

        READ_PROCESSING(
            ch_inputs,
            ch_align_dna_tumor_out,
            ch_inputs.map { meta -> [meta, [], []] }, // ch_normal_bam
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.unmap_regions,
            false,  // has_umis
        )

        ch_versions = ch_versions.mix(READ_PROCESSING.out.versions)

        ch_process_dna_tumor_out = ch_process_dna_tumor_out.mix(READ_PROCESSING.out.dna_tumor)

    } else {

        ch_process_dna_tumor_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = Channel.empty()
    if (run_config.stages.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ch_process_dna_tumor_out,
            ch_process_dna_normal_out,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
            [],  // target_region_bed
        )

        ch_versions = ch_versions.mix(AMBER_PROFILING.out.versions)

        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)

    } else {

        ch_amber_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = Channel.empty()
    if (run_config.stages.cobalt) {

        COBALT_PROFILING(
            ch_inputs,
            ch_process_dna_tumor_out,
            ch_inputs.map { meta -> [meta, [], []] }, // ch_normal_bam
            hmf_data.gc_profile,
            hmf_data.diploid_bed,
            [],  // panel_target_region_normalisation
        )

        ch_versions = ch_versions.mix(COBALT_PROFILING.out.versions)

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, []] }

    }






    // TODO(SW): rework to accept multiple samples per grouping (patient/subject)

    //
    // SUBWORKFLOW: Append new sample data to primary SAGE WGS VCF
    //
    // channel: [ meta, sage_append_vcf ]
    ch_sage_somatic_append_out = Channel.empty()
    if (run_config.stages.orange) {

        SAGE_APPEND(
            ch_inputs,
            ch_inputs.map { meta -> [meta, []] },  // ch_purple_dir
            ch_process_dna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
        )

        ch_versions = ch_versions.mix(SAGE_APPEND.out.versions)
        ch_sage_somatic_append_out = ch_sage_somatic_append_out.mix(SAGE_APPEND.out.somatic_vcf)

    } else {

        ch_sage_somatic_append_out = ch_inputs.map { meta -> [meta, []] }

    }






    // TODO(SW): complete subworkflow implementation

    //
    // SUBWORKFLOW: Run WISP to estimate tumor purity
    //
    if (run_config.stages.wisp) {

        WISP_ANALYSIS(
            ch_inputs,
            ch_amber_out,
            ch_cobalt_out,
            ch_sage_somatic_append_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
        )

        ch_versions = ch_versions.mix(WISP_ANALYSIS.out.versions)

    }

    //
    // TASK: Aggregate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
