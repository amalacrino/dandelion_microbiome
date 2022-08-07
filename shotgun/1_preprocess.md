### 1. TrimGalore

```bash
$DATADIR=add_path_here
$OUTDIR=add_path_here
find -name "*_R1.fq.gz" | cut -d "_" -f1 | parallel -j 16 trim_galore --illumina --paired --fastqc --gzip -o $OUTDIR/ {}\_R1.fq.gz {}\_R2.fq.gz
```
