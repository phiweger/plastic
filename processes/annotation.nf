process map_proteins {
    container 'nanozoo/miniprot:2.24--6d307d9'
    cpus 8
    debug true

    input:
        tuple(val(name), path(genome), path(proteins))

    output:
        tuple(val(name), path(genome), path('pet.paf'))

    """
    miniprot --outn 1 -S -t ${task.cpus} ${genome} ${proteins} > pet.paf
    """
}