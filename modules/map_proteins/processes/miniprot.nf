process miniprot {
    container 'nanozoo/miniprot:2.24--6d307d9'
    cpus 1
    // debug true

    input:
        tuple(val(name), path(genome), path(proteins))

    output:
        tuple(val(name), path(genome), path("${name}.paf"))

    """
    miniprot --outn 1 -S -t ${task.cpus} ${genome} ${proteins} > ${name}.paf
    """
}
