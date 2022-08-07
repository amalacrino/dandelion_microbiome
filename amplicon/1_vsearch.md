Data handling was carried out using VSEARCH (Rognes et al. 2016), assigning taxonomy using SILVA v132 (16S; Quast et al. 2012) or UNITE () reference databases.

Run the pipeline according to https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline and assign taxonomy using the command

```bash
vsearch --sintax $PRJ_DIR/otus.fasta --db $PRJ_DIR/silva_16s_v131.fa --tabbedout $PRJ_DIR/ASV_tax_raw.txt --sintax_cutoff 0.5 --threads $NTHREADS
```
