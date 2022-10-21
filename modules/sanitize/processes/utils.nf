process rename {
    cache 'lenient'  // TODO: not sure why this is needed, file sys changes?
    container 'nanozoo/sourmash:4.2.3--ce062ab'
    // publishDir "${params.results}/names", pattern: '*.json', mode: 'copy'
    // errorStrategy 'ignore'
    // debug true
    /*
    In strict mode the process will fail if a duplicate contig is detected,
    otherwise the contig is skipped and thus not considered in further
    analyses.

    rename.py will automatically read in gzipped and unzipped files, and
    write out uncompressed ones, because there are downstream programs like
    "prodigal" and "prokka" which cannot process compressed ones. So this is
    unnecessary:

    if [[ ${genome} =~ \\.gz\$ ]]
        then
            gunzip -c ${genome} > tmp
        else
            cat ${genome} > tmp
    fi
    */

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path("${name}.fna.gz"), path("${name}.json"))

    script:
        if( params.strict )
            """
            rename.py --name ${name} --genome ${genome} --strict
            """

        else
            """
            rename.py --name ${name} --genome ${genome}
            gzip ${name}.fna
            """

}


process checksum {
    cache 'lenient'  // path and size-based
    // cache 'deep'  // content-based but slow

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), env(hash), path(genome))

    shell:
    '''
    hash=$(md5sum !{genome} | awk '{ printf $1 }')
    '''
    // awk print will insert newline, printf does not
    // stackoverflow.com/questions/2021982
    // alternative: ... | tr -d '\n'
    // stackoverflow.com/questions/3134791
}


process merge_names {
    publishDir "${params.results}", mode: 'copy'
    
    input:
        path(names)

    output:
        path('names.json')

    """
    jq -s '.' ${names} > names.json
    """
}


process minlen {
    container 'nanozoo/seqkit:2.3.1--401ce9e'
    cpus 8

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path("${name}.filt.fna.gz"))

    """
    seqkit seq \
        --seq-type dna \
        --line-width 80 \
        --upper-case \
        --remove-gaps \
        --min-len ${params.minlen} \
        --validate-seq \
        --threads ${task.cpus} \
        ${genome} \
        > tmp
    gzip -c tmp > ${name}.filt.fna.gz
    """
}
