nextflow.enable.dsl = 2
// include {sayHello} from './some/module' addParams(foo: 'Ciao')
// https://www.nextflow.io/docs/latest/dsl2.html#module-aliases


/*
nextflow run plastic/main.nf --genomes gtdb_genomes_reps_r207 --proteins plastic/data/pet.faa --db gtdb-rs207.genomic.k21.lca.json.gz --results results -resume --debug
*/


include { 
    sanitize as clean_genomes
    sanitize as clean_proteins
    } from './modules/sanitize/main.nf'


include { taxa } from './modules/taxa/main.nf'
include { map_proteins } from './modules/map_proteins/main.nf'


workflow {
    
    // Oftentimes, genome files have associated eg sample names; so we proceed
    // with the (name, path) pattern here.
    suffix = '*{.fna,.fna.gz,.fasta,.fasta.gz}'
    genomes = channel.fromPath("${params.genomes}/${suffix}")
                     .map { x -> [x.baseName, x] }

    db = channel.fromPath(params.db)

    proteins = channel.fromPath("${params.proteins}")
    /*
    TODO: Clean protein file, too

    foo = channel.fromPath("${params.proteins}")
                      .map { x -> [x.baseName, x] }
    clean_proteins(foo)
    */

    // TODO: rename, but use the https://github.com/nextflow-io/nf-sqldb
    if (params.debug)
        genomes.take(200) | clean_genomes
        // .randomSample() is not chached
    else
        genomes | clean_genomes

    map_proteins(clean_genomes.out.sequences, proteins)

    taxa(map_proteins.out, db)

}


