include { pick_best; save_pick } from './processes/pick.nf'
include { miniprot } from './processes/miniprot.nf'


workflow map_proteins {
    take:
        genomes
        proteins

    main:
        best = genomes | combine(proteins) | miniprot | pick_best
        groups = best | groupTuple | map { i -> [i[0], i[4], i[5]] } | transpose
        groups | save_pick

    emit:
        pick_best.out.map { i -> [i[1], i[2]] }
}