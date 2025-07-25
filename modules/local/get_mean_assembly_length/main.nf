process GET_MEAN_ASSEMBLY_LENGTH {

    label 'process_single'

    tag "$family"

    errorStrategy = {
        if (task.exitStatus == 100) {
            // ignoring cases when family does not have children
            log.warn("Could not find species for family ${family}.")
            return 'ignore'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/7140956576e6779ff2bea610b8e41bde95d4909d67b31634d0ddb6dba50aef5a/data':
        'community.wave.seqera.io/library/requests_tenacity_xmltodict:9e74a2aeeb88aab9' }"

    input:
    val family

    output:
    tuple val(family), eval('cat mean_assembly_length.txt'),                                                                     emit: families
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'),            topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_assembly_stats.py --family $family
    """

    stub:
    """
    echo 1 > mean_assembly_length.txt
    """

}
