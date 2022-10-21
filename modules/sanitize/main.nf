include { rename; checksum; merge_names; minlen } from './processes/utils.nf'


workflow sanitize {
    take:
        genomes

    main:
        genomes | checksum | map { i -> [i[1], i[2]] } | rename

        rename.out | map { i -> i[2] } | collect | merge_names
        rename.out | map { i -> [i[0], i[1]] } | minlen

    emit:
        // https://www.nextflow.io/docs/edge/dsl2.html#workflow-named-output
        sequences = minlen.out
        names = merge_names.out
}


// TODO: We should be able to test this module standalone
workflow {
    genomes = channel.fromPath('test/*{.fna,.fna.gz}')
                     .map { x -> [x.baseName, x] }

    genomes | sanitize

}