include { 
    rename
    checksum
    } from '../processes/utils.nf'


workflow wf_sanitize {
    take:
        genomes

    main:
        cs = checksum(genomes).map { it -> [it[1], it[2]] }
        rename(cs)
    emit:
        // https://www.nextflow.io/docs/edge/dsl2.html#workflow-named-output
        cleanomes = rename.out.map(it -> [it[0], it[1]])
        name_map = rename.out.map(it -> [it[2]])
}