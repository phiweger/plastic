include { sketch; classify } from './processes/sourmash.nf'


workflow taxa {
    take:
        genomes
        db

    main:
        sigs = genomes | sketch | map { i -> i[2] } | collect
        classify(sigs, db)

    emit:
        classify.out

}