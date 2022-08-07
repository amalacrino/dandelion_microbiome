Data handling was carried out using VSEARCH v2.21.1 (Rognes et al. [2016](https://peerj.com/articles/2584/)), assigning taxonomy using SILVA v132 (16S; Quast et al. [2012](https://academic.oup.com/nar/article/41/D1/D590/1069277)) or UNITE v8.3 (ITS; Nilsson et al. [2019](https://academic.oup.com/nar/article/47/D1/D259/5146189)) reference databases.

We ran the pipeline according to https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline, and assigned taxonomy using the command:

```bash
vsearch --sintax $PRJ_DIR/otus.fasta --db $PRJ_DIR/silva_16s_v131.fa --tabbedout $PRJ_DIR/ASV_tax_raw.txt --sintax_cutoff 0.5 --threads $NTHREADS
```
