nextflow.enable.dsl = 2
// include {sayHello} from './some/module' addParams(foo: 'Ciao')
// https://www.nextflow.io/docs/latest/dsl2.html#module-aliases


/*
nextflow run plastic/main.nf --genomes gtdb_genomes_reps_r207 --proteins plastic/data/pet.faa --db gtdb-rs207.genomic.k21.lca.json.gz --results results -resume --debug
*/

include { 
    map_proteins
    } from './processes/annotation.nf'


include { 
    sanitize as clean_genomes
    sanitize as clean_proteins
    } from './modules/sanitize/main.nf'


process filt {
    // debug true
    errorStrategy 'ignore'
    container 'nanozoo/miniprot:2.24--0c673d2'
    publishDir "${params.results}/proteins", pattern: '*.faa', mode: 'copy'
    publishDir "${params.results}/coding", pattern: '*.fna', mode: 'copy'

    input:
        tuple(val(name), path(genome), path(aln))

    output:
        tuple(
            env(QRY),
            val(name),
            path(genome),
            path(aln),
            path("${name}.best_hit.fna"),
            path("${name}.best_hit.faa"),
            ) optional true

    script:
    """
    QRY=\$(filter_protein_aln.py --identity 0.7 --coverage 0.7 --aln ${aln} --genome ${genome} --out ${name}.best_hit)
    """
}


process sketch {
    // debug true
    container 'nanozoo/sourmash:4.2.3--ce062ab'
    publishDir "${params.results}/genomes", pattern: '*.fna.gz', mode: 'copy'
    publishDir "${params.results}/aln", pattern: '*.paf', mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(group), val(name), path(genome), path(aln), path(fna), path(faa))

    output:
        tuple(val(name), path(genome), path(aln), path("${name}.sig"))

    """
    sourmash sketch dna -p 'k=21,scaled=1000' ${genome} -o ${name}.sig
    """
}


process taxa {
    container 'nanozoo/sourmash:4.2.3--ce062ab'
    publishDir "${params.results}", mode: 'copy'

    input:
        path(sketches)
        path(db)

    output:
        path('taxa.csv')

    """
    sourmash lca classify --db ${db} --query ${sketches} -o taxa.csv
    """
}


workflow {
    
    // Oftentimes, genome files have associated eg sample names; so we proceed
    // with the (name, path) pattern here.
    genomes = channel.fromPath("${params.genomes}/*{.fna,.fna.gz}")
                     .map { x -> [x.baseName, x] }
    // genomes = channel.fromPath("${params.genomes}/**/*{.fna,.fna.gz}")
    db = channel.fromPath(params.db)

    proteins = channel.fromPath("${params.proteins}")

    /*
    TODO: Clean protein file, too

    foo = channel.fromPath("${params.proteins}")
                      .map { x -> [x.baseName, x] }
    clean_proteins(foo)
    */

    if (params.debug)
        genomes.take(200) | clean_genomes
        // .randomSample() is not chached
    else
        genomes | clean_genomes

    
    // name_map = wf_sanitize.out.name_map    
    // proteins | wf_sanitize

    clean_genomes.out.sequences | combine(proteins) | map_proteins | filt

    //filt.out.view()
    
    //filt.out.view()
    //sketch(filt.out)
    //filt(map_proteins.out)

    // fn = { it -> [it[0], it[1..-1]]}
    // x = filt.out.map(fn)
    // y = filt.out
    // x.join(x)
    // y.join(y, by: 0).view()
    filt.out.groupTuple().map { it -> [it[0], it[1]] }.transpose().view()
    // transpose() i opposite of groupTuple
    // filt.out.join().view()
    
    // .map{it -> [it[1], it[2]]}

    sketch(filt.out)
    // taxa(sketch.out.map { it -> [it[3]] }.collect(), db)

}


