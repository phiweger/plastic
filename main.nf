nextflow.enable.dsl = 2

/*
nextflow run plastic/main.nf --genomes gtdb_genomes_reps_r207 --proteins plastic/data/pet.faa --db gtdb-rs207.genomic.k21.lca.json.gz --results results -resume --debug
*/

include { 
    map_proteins
    } from './processes/annotation.nf'


include { 
    wf_sanitize
    } from './workflows/sanitize.nf'


process filt {
    // debug true
    errorStrategy 'ignore'
    container 'nanozoo/miniprot:2.24--6d307d9'

    input:
        tuple(val(name), path(genome), path(aln))

    output:
        tuple(val(name), path(genome), path(aln), env(FOUND))

    script:
    """
    FOUND=\$(filter_protein_aln.py --identity 0.8 --aln ${aln})
    """
}


process sketch {
    // debug true
    container 'nanozoo/sourmash:4.2.3--ce062ab'
    publishDir "${params.results}/targets", pattern: '*.fna.gz', mode: 'copy'
    publishDir "${params.results}/aln", pattern: '*.paf', mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(name), path(genome), path(aln), val(found))

    output:
        tuple(val(name), path(genome), path(aln), path("${name}.sig"))

    script:
        if ( found == '1' )
            """
            # echo ${name}
            sourmash sketch dna -p 'k=21,scaled=1000' ${genome} -o ${name}.sig
            """
        // else
        //     // """
        //     // echo ${genome}
        //     // """
        //     error "Not found."
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
    genomes = channel.fromPath("${params.genomes}/**/*{.fna,.fna.gz}")
                     .map { x -> [x.baseName, x] }

    db = channel.fromPath(params.db)

    proteins = channel.fromPath("${params.proteins}")
    
    if (params.debug)
        genomes.take(5000) | wf_sanitize
        // .randomSample() is not chached
    else
        genomes | wf_sanitize

    
    // name_map = wf_sanitize.out.name_map    

    wf_sanitize.out.cleanomes.combine(proteins) | map_proteins | filt | sketch
    //filt.out.view()
    //sketch(filt.out)
    //filt(map_proteins.out)

    //sketch(filt.out)
    taxa(sketch.out.map { it -> [it[3]] }.collect(), db)

}


