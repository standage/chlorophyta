```bash
genhub-build.py --cfgdir=cfg --batch=chlorophyta+ \
                --workdir=species --numprocs=16 \
                download format prepare stats cluster

python status.py GenHub.hiloci.tsv > Chlorophyta.hiLocus.status.tsv

python breakdown.py <(cat species/*/*.iloci.tsv) Chlorophyta.hiLocus.status.tsv \
    > Chlorophyta-breakdown-bp.tsv
python breakdown.py --counts <(cat species/*/*.iloci.tsv) \
                    Chlorophyta.hiLocus.status.tsv \
    > Chlorophyta-breakdown-counts.tsv
```
