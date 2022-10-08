## README

Find all genomes that carry a copy of a set of target proteins. Use case here: 
Find PETases in GTDB genomes.

```bash
nextflow run main.nf \
    --genomes genomes \
    --proteins pet.faa \
    --db db/gtdb-rs207.genomic.k21.lca.json.gz \
    --results results \
    -resume
```
