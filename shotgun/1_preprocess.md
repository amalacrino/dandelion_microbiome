### 1. TrimGalore

```bash
DATADIR=$add_path_here
OUTDIR=$add_path_here
nthreads = 16

cd $DATADIR
find -name "*_R1.fq.gz" | cut -d "_" -f1 | parallel -j $nthreads trim_galore --illumina --paired --fastqc --gzip -o $OUTDIR/ {}\_R1.fq.gz {}\_R2.fq.gz
```

### 2. Remove host plant reads

```bash
DATADIR=$add_path_here
REFDIR=$add_path_here
OUTDIR=$add_path_here
nthreads = 16

cd $REFDIR
wget ftp://download.big.ac.cn/gwh/Plants/Taraxacum_kok-saghyz_TKS_GWHAAAA00000000/GWHAAAA00000000.genome.fasta.gz
gunzip GWHAAAA00000000.genome.fasta.gz
bowtie2-build --threads 8 GWHAAAA00000000.genome.fasta dandelion

cd $DATADIR
for file in *_R1_val_1.fq.gz
do
    SAMPLE=`basename $file _R1_val_1.fq.gz`
    echo $SAMPLE Bowtie2
    bowtie2 -p 32 -x $REFDIR/dandelion -1 ${SAMPLE}_R1_val_1.fq.gz -2 ${SAMPLE}_R2_val_2.fq.gz -S $OUTDIR/${SAMPLE}_mapped_and_unmapped.sam
    echo $SAMPLE conversion from sam to bam
    samtools view -@ 32 -bS $OUTDIR/${SAMPLE}_mapped_and_unmapped.sam > $OUTDIR/${SAMPLE}_mapped_and_unmapped.bam
    echo $SAMPLE select unmapped
    samtools view -@ 32 -b -f 12 -F 256 $OUTDIR/${SAMPLE}_mapped_and_unmapped.bam > $OUTDIR/${SAMPLE}_bothEndsUnmapped.bam
    echo $SAMPLE sorting unmapped
    samtools sort -@ 32 -n $OUTDIR/${SAMPLE}_bothEndsUnmapped.bam -o $OUTDIR/${SAMPLE}_bothReadsUnmapped_sorted.bam
    echo $SAMPLE extract fastq
    samtools fastq -@ 32 $OUTDIR/${SAMPLE}_bothReadsUnmapped_sorted.bam -1 $OUTDIR/${SAMPLE}_host_removed_R1.fastq.gz -2 $OUTDIR/${SAMPLE}_host_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n
    find $OUTDIR -type f  ! -name "*.fastq.gz"  -delete
done

```
