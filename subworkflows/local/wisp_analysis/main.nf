//
// WISP estimates tumor purity in longitudinal samples using WGS data of the primary
//

import Constants
import Utils

include { WISP } from '../../../modules/local/wisp/main'

workflow WISP_ANALYSIS {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_amber_out               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_out              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_sage_somatic_append_out // channel: [mandatory] [ meta, sage_append_vcf ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_fai       // channel: [mandatory] /path/to/genome_fai

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, ... ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ...,
    )
        .map { ... ->
            return [ ... ]
        }
        .branch { ... ->
            runnable: ...
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_wisp, ... ]
    ch_amber_inputs = ch_inputs_sorted.runnable
        .map { meta, ... ->
            return [meta, ...]
        }

    // Run process
    WISP(
    )

    ch_versions = ch_versions.mix(AMBER.out.versions)

    emit:
    versions     = ch_versions     // channel: [ versions.yml ]
}
