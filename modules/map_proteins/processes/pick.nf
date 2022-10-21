process pick_best {
    // debug true
    errorStrategy 'ignore'
    container 'nanozoo/miniprot:2.24--0c673d2'
    // publishDir "${params.results}/proteins/${match}", pattern: '*.faa', mode: 'copy'
    // publishDir "${params.results}/coding/${match}", pattern: '*.fna', mode: 'copy'

    input:
        tuple(val(name), path(genome), path(aln))

    output:
        tuple(
            env(match),
            val(name),
            path(genome),
            path(aln),
            path("${name}.best_hit.fna"),
            path("${name}.best_hit.faa"),
            ) optional true

    script:
    """
    match=\$(filter_protein_aln.py --identity 0.7 --coverage 0.7 --aln ${aln} --genome ${genome} --out ${name}.best_hit)
    """
}


process save_pick {
    /*
    This is rather dumb, but the match env variable from the "filt" process
    is not available for the publishDir directive, and so this process exists
    with the only task of organising the results.

    Related discussion, but AFAIK unresolved:

    https://github.com/nextflow-io/nextflow/discussions/1933
    */
    publishDir "${params.results}/proteins/${match}", pattern: '*.faa', mode: 'copy'
    publishDir "${params.results}/coding/${match}", pattern: '*.fna', mode: 'copy'

    input:
        tuple(val(match), path(gene), path(protein))

    output:
        tuple(val(match), path(gene), path(protein))

    """
    sleep 1
    """
}
