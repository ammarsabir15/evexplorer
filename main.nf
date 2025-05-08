#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/evexplorer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/evexplorer
    Website: https://nf-co.re/evexplorer
    Slack  : https://nfcore.slack.com/channels/evexplorer
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EVEXPLORER  } from './workflows/evexplorer'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_evexplorer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_evexplorer_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_evexplorer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
     params.fasta = getGenomeAttribute('fasta')
//     ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
//                         .map { file -> tuple('genome', file) }
//


params.chr_names = '/references/DERFINDER_ref/STAR_INDEX/chrName.txt'
params.gtf_1 = '/references/Derfinder_pipeline/gencodeV38_pluspirnadb1_7_6_plusmirbase21_tRNAscan_MINT_hg38_V1.gff3'
params.gtf_2 = '/references/Derfinder_pipeline/homo_sapiens.GRCh38.gff3'

// Define channels for GTF files
ch_gtf_1 = Channel.fromPath(params.gtf_1)
ch_gtf_2 = Channel.fromPath(params.gtf_2)

// Define channel for chromosome names
ch_chr_names = Channel.fromPath(params.chr_names)


ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
                        .map { file -> tuple('genome', file) }
                        .first()
                        .map { id, file -> tuple(id, file) }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_EVEXPLORER {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    EVEXPLORER (
        samplesheet,
        ch_fasta,
        ch_gtf_1,
        ch_gtf_2,
        ch_chr_names
)

    emit:
    multiqc_report = EVEXPLORER.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.resume
 )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_EVEXPLORER (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_EVEXPLORER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


