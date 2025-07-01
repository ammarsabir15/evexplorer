#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    nf-core/evexplorer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Github : https://github.com/nf-core/evexplorer

    Website: https://nf-co.re/evexplorer
    Slack  : https://nfcore.slack.com/channels/evexplorer
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
*/


include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_evexplorer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_evexplorer_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_evexplorer_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
*/


params.fasta = getGenomeAttribute('fasta')
params.chr_names = '/references/DERFINDER_ref/STAR_INDEX/chrName.txt'
params.gtf_1 = '/references/Derfinder_pipeline/gencodeV38_pluspirnadb1_7_6_plusmirbase21_tRNAscan_MINT_hg38_V1.gff3'
params.gtf_2 = '/references/Derfinder_pipeline/homo_sapiens.GRCh38.gff3'

// Channels
ch_gtf_1     = Channel.fromPath(params.gtf_1)
ch_gtf_2     = Channel.fromPath(params.gtf_2)
ch_chr_names = Channel.fromPath(params.chr_names)
ch_fasta     = Channel.fromPath(params.fasta, checkIfExists: true)
                  .map { file -> tuple('genome', file) }
                  .first()
                  .map { id, file -> tuple(id, file) }


    include { SMALL_READS } from './workflows/small_reads'
    include { LONG_READS } from './workflows/long_reads'

	
//
// WORKFLOW: Run main analysis pipeline depending on type of input
//


workflow NFCORE_EVEXPLORER {
    take:
    samplesheet
	
	
    main:
	//
    // WORKFLOW: Run pipeline
    //
    if (params.platform == 'comboseq' || params.platform == 'nextflex') {
        SMALL_READS (
	samplesheet,
        ch_fasta,
        ch_gtf_1,
        ch_gtf_2,
        ch_chr_names
		)
		multiqc_report_ch = SMALL_READS.out.multiqc_report
    } else if (params.platform == 'ont') {
        LONG_READS (
        samplesheet,
        ch_fasta,
        ch_gtf_1,
        ch_gtf_2,
        ch_chr_names		
		)
		multiqc_report_ch = LONG_READS.out.multiqc_report
    }
	emit:
    multiqc_report = multiqc_report_ch
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//

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
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_EVEXPLORER (
        PIPELINE_INITIALISATION.out.samplesheet,
	)
	
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
