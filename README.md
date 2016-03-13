# Analysis of 9 chlorophyte genomes

## Data retrieval, iLocus parsing, and protein clustering

We use GenHub to download the 9 chlorophyte genomes (as well as 4 outgroups) from RefSeq, compute iLoci, and cluster iLocus protein products from all species.
We then use a custom Python script to assign each iLocus a provisional status based on the proteins with which it clusters.
- *HighlyConserved*: iLocus has a protein product conserved in all 9 chlorophyte species
- *Conserved*: iLocus has a protein product conserved in at least 4 of the 9 chlorophyte species
- *Matched*: iLocus has a protein product that is conserved in at least on other species, including the outgroups
- *Unmatched*: iLocus has a protein product that is not conserved in any other species

```bash
genhub-build.py --cfgdir=cfg --batch=chlorophyta+ \
                --workdir=species --numprocs=16 \
                download format prepare stats cluster
python status.py GenHub.hiloci.tsv > Chlorophyta.hiLocus.pre-status.tsv
```

## Post-processing of unmatched iLoci

The clustering criteria place somewhat strict limits on length differences between proteins.
Many iLoci labeled *Unmatched* by the initial step are therefore not truly orphans, but may be annotated incompletely or incorrectly in one or more species, or may have evolved to such an extent that near-full-length alignments are not possible.
Most of these *Unmatched* iLoci will indeed encode *bona fide* protein products that are likely conserved in other species, but may need additional attention before they can be reliably used for comparative genomics analysis.
This post-processing procedure will distinguish these initially unmatched iLoci from iLoci that have no reliable matches in other species, the former being relabeled as *Matched* and the latter being relabeled as *Orphan*.

```bash
# Build BLAST search index
mkdir blastdbs
for species in Apro Crei Csub Cvar Mpus Msrc Oluc Otau Vcar 
do
    makeblastdb -in species/${species}/${species}.prot.fa -dbtype prot -parse_seqids
done

for species in Apro Crei Csub Cvar Mpus Msrc Oluc Otau Vcar
do
    # Grab proteins from iLoci initially labeled  "Unmatched"
    grep Unmatched Chlorophyta.hiLocus.pre-status.tsv \
        | grep $species \
        | cut -f 5 \
        > blastdbs/${species}.unmatched.txt
    blastdbcmd -db species/${species}/${species}.prot.fa \
               -entry_batch blastdbs/${species}.unmatched.txt \
        > blastdbs/${species}.unmatched.fa

    # Compose a database of proteins from all species except this one
    blastdb_aliastool -title ${species} \
                      -out blastdbs/${spec} \
                      -dblist_file <(ls species/*/????.prot.fa | grep -v $spec) \
                      -dbtype prot

    # Execute the BLAST search
    blastp -query blastdbs/${species}.unmatched.fa \
           -db blastdbs/${species} \
           -num_alignments 5 \
           -num_descriptions 5 \
           -evalue 1e-4 \
           -out blastdbs/${species}.unmatched.blastp

    # Extract protein IDs of newly matched proteins
    MuSeqBox -i blastdbs/${species}.unmatched.blastp \
             -L 100 -d 16 -l 24 -c crtfile \
        > blastdbs/${species}.unmatched.msb
    grep -e '^XP' -e '^NP' blastdbs/${species}.unmatched.msb \
        | awk '{ print $1 }' \
        | sort \
        | uniq \
        > blastdbs/${species}.blast_matches.txt
done

# Re-classify "Unmatched" iLoci as "Matched" or "Orphan"
python post_blast.py <(cat blastdbs/*.blast_matches.txt) Chlorophyta.hiLocus.pre-status.tsv \
    > Chlorophyta.hiLocus.status.tsv
```

Finally, we use custom Python scripts to compute a breakdown of each species, for each species showing the number of iLoci assigned to each category, as well as the proportion of the genome occupied by iLoci of that category.
See [notebook.ipynb](notebook.ipybn) for the code used to plot these breakdowns.

```bash
python breakdown.py <(cat species/*/*.iloci.tsv) Chlorophyta.hiLocus.status.tsv \
    > Chlorophyta-breakdown-bp.tsv
python breakdown.py --counts <(cat species/*/*.iloci.tsv) \
                    Chlorophyta.hiLocus.status.tsv \
    > Chlorophyta-breakdown-counts.tsv
```
