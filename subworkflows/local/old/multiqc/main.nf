include { MULTIQC                                                 } from '../../../modules/nf-core/multiqc'

include { paramsSummaryMap                                        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                    } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                  } from '../../nf-core/utils_nfcore_pipeline'
include { formatVersionsToYAML                                    } from '../utils_nfcore_nf-chimera-detector_pipeline'
include { methodsDescriptionText                                  } from '../utils_nfcore_nf-chimera-detector_pipeline'


workflow MULTIQC_WORKFLOW {

    take:
    ch_multiqc_files
    ch_versions

    main:


    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    // Collate and save software versions obtained from topic channels
    // TODO: use the nf-core functions when they are adapted to channel topics

    // Collate and save software versions
    formatVersionsToYAML ( Channel.topic('versions') )
        .mix ( softwareVersionsToYAML( ch_versions ) ) // mix with versions obtained from emit outputs
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    // ------------------------------------------------------------------------------------
    // CONFIG
    // ------------------------------------------------------------------------------------

    ch_multiqc_config = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )

    summary_params = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value( paramsSummaryMultiqc(summary_params) )

    ch_multiqc_custom_config = params.multiqc_config ?
                                    Channel.fromPath(params.multiqc_config, checkIfExists: true) :
                                    Channel.empty()

    ch_multiqc_logo = params.multiqc_logo ?
                        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
                        Channel.empty()

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
                                                file(params.multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description = Channel.value( methodsDescriptionText(ch_multiqc_custom_methods_description) )

    // Adding metadata to MultiQC
    ch_multiqc_files = ch_multiqc_files
                            .mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
                            .mix( ch_collated_versions )
                            .mix( ch_methods_description.collectFile( name: 'methods_description_mqc.yaml', sort: true ) )

    // ------------------------------------------------------------------------------------
    // DATA
    // ------------------------------------------------------------------------------------

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )


    emit:
    multiqc_report = MULTIQC.out.report
}

