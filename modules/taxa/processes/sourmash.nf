process sketch {
    // debug true
    container 'nanozoo/sourmash:4.2.3--ce062ab'
    publishDir "${params.results}/genomes", pattern: '*.fna.gz', mode: 'copy'
    // publishDir "${params.results}/aln", pattern: '*.paf', mode: 'copy'
    // errorStrategy 'ignore'

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path(genome), path("${name}.sig"))

    """
    sourmash sketch dna -p 'k=21,scaled=1000' ${genome} -o ${name}.sig
    """
}


process classify {
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
