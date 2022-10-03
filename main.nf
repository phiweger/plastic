nextflow.enable.dsl = 2


include { 
    map_proteins
    } from './processes/annotation.nf'


include { 
    wf_sanitize
    } from './workflows/sanitize.nf'


process filt {
    debug true
    container 'nanozoo/miniprot:2.24--6d307d9'

    input:
        tuple(val(name), path(genome), path(aln))

    output:
        tuple(val(name), path(genome), env(FOUND))

    script:
    """
    FOUND=\$(filter_protein_aln.py --identity 0.8 --aln ${aln})
    """
}


process sketch {
    debug true
    container 'nanozoo/sourmash:4.2.3--ce062ab'
    publishDir "${params.results}/targets/*.fna", mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(name), path(genome), val(found))

    output:
        tuple(val(name), path(genome), path("${name}.sig"))

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
    genomes = channel.fromPath("${params.genomes}/*{.fna,.fna.gz}")
                     .map { x -> [x.baseName, x] }

    db = channel.fromPath(params.db)

    proteins = channel.fromPath("${params.proteins}")
    
    if (params.debug) {
        checksum(genomes.take(20))
        // .randomSample() is not chached
    }

    genomes | wf_sanitize
    // name_map = wf_sanitize.out.name_map    

    wf_sanitize.out.cleanomes.combine(proteins) | map_proteins | filt | sketch
    //filt(map_proteins.out)

    //sketch(filt.out)
    taxa(sketch.out.map { it -> [it[2]] }.collect(), db)

}


