### TrimGalore

We used TrimGalore to clean raw data from adaptors and low quality reads.

```bash
DATADIR=$add_path_here
OUTDIR=$add_path_here
nthreads = 16

cd $DATADIR
find -name "*_R1.fq.gz" | cut -d "_" -f1 | parallel -j $nthreads trim_galore --illumina --paired --fastqc --gzip -o $OUTDIR/ {}\_R1.fq.gz {}\_R2.fq.gz
```

### Remove host plant reads

We used bowtie to map raw data to the *T. koksaghyz* reference genome. The, we used samtools to select reads mapping the host plant genome and convert them back into fastq files to be used for downstream analyses.

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

### MegaHit

We assembled the raw reads into contigs using MegaHit. From this experiment we obtained a massive amount of data, and MegaHit requires too much memory to process it all together. So, here, we are going to perform the contig assebly separately for each sample and later in the pipeline we will merge this output. Using the raw reads filtered from the host reads as input, we can run this code:

```bash
DATADIR=$add_path_here
OUTDIR=$add_path_here

cd $DATADIR

for file in *_host_removed_R1.fastq.gz
do
    SAMPLE=`basename $file _host_removed_R1.fastq.gz`
		megahit -t 32 -1 ${SAMPLE}_host_removed_R1.fastq.gz \
			-2 ${SAMPLE}_host_removed_R2.fastq.gz \
			-o $OUTDIR/${SAMPLE}
done 
```

MegaHit creates a folder for each of our inputs, with a file named `final.contigs.fa` were it stores the information we will need for the nexrt steps. So, running the code below we will take each `final.contigs.fa` file, rename it according to the folder name (which is our sample), and move it to the new folder.

```bash
MegahitPATH=$add_path_here
MegahitFILES=$add_path_here
find "$MegahitPATH" -type f -iname 'final.contigs.fa' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" "$MegahitFILES/${path##*/}.fa" ' sh_cp {} \;

cd $MegahitFILES
cat *.fa > allseqs.fa
```

MegaHit gives the same names to the contigs obtained from different files (something like contig1. contig2, ...) and some softwares down in the pipeline do not like to have contigs with the same name within the same `fasta` file. So, here we are going to rename all the contigs in the file we just obtained. First, we need to know how many contigs our file has:

```bash
sed '/^>/d' allseqs.fa | wc -l
```
In our case the result was `29853253`. So we are going to use this number to run the line below.

```bash
for i in {1..29853253}; do echo A_$i; done | paste - <(sed '/^>/d' allseqs.fa) | sed -e 's/^/>/' -e 's/\t/\n/' > new_file.fa
```

### Remove duplicate reads

We removed duplicate contigs to avoid any interference with the pipeline below.

```bash
seqkit rmdup -s < new_file.fa > input.faa

rm allseqs.fa new_file.fa
```

### Prokka

And we finally ran Prokka for functional annotation of the contigs identified above. Given that Prokka was not working well with a single file with all the contigs, we optimized this process by breaking it in several files with 10,000 contigs each and processing them in parallel.

```bash
DATADIR=$add_path_here
OUTDIR=$add_path_here

cd $DATADIR
cp input.faa $OUTDIR
cd $OUTDIR

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < input.faa

find -name "*.fa" | parallel -j 32 prokka --outdir $OUTDIR --norrna --notrna --metagenome --force --cpus 1 --prefix {}\.prokka {}
```

Then we collated all `*.tsv` files all together:

```bash
DATADIR=$add_path_here
OUTDIR=$add_path_here

cd $DATADIR
cp *.tsv $OUTDIR
cd $OUTDIR
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.tsv > prokka.txt
rm *.tsv
```

### Counts table

Given that we wanted to test the effect of our treatments on the frequency of genes within the microbiome, we used Bowtie to map each sample to the contig reference, generating a table with counts for each contig and sample.


```bash
DATADIR=$add_path_here
RAWREADS=$add_path_here
INDIR=$add_path_here
OUTDIR=$add_path_here
OUTDIR2=$add_path_here

cd $DATADIR
cat *.ffn > prokka.fasta
cp prokka.fasta $INDIR

cd $INDIR

bowtie2-build --threads 16 prokka.fasta reference
samtools faidx prokka.fasta

cd $RAWREADS

for file in *_R1_val_1.fq.gz
do
  FILENAME=`basename $file _R1_val_1.fq.gz`;
  echo $FILENAME;
  bowtie2 -x $INDIR/reference -1 ${FILENAME}_R1_val_1.fq.gz -2 ${FILENAME}_R2_val_2.fq.gz -p 16 -S $OUTDIR/${FILENAME}.sam;
done

cd $OUTDIR

cp *.sam $OUTDIR2

cd $OUTDIR2

for file in *.sam
do
  f=`basename $file .sam`;
  echo $f;
  samtools view -@ 16 -bS ${f}.sam > ${f}.bam;
  samtools sort -@ 16 ${f}.bam -o ${f}.sorted;
  samtools index -@ 16 ${f}.sorted;
  samtools idxstats ${f}.sorted > ${f}.idxstats.txt;
done

```

And to generate the count table we used this script in Python 2.7 (`get_count_table.py`) obtained from [here](https://github.com/edamame-course/Metagenome/blob/master/get_count_table.py)

```bash
#!/usr/bin/python
#this script merge idxstat file into one
#usage: python get_count_table.py *.idxstats.txt > count.txt
#result: contig_id, length, counts

import sys
flag = 0
ids = []
length = []
count = []
fname = []
for f in sys.argv[1:]:
    fname.append(f)
    co = []
    for n,line in enumerate(open(f)):
        #print line
        if(line[:1] == '*'):
            continue
        spl = line.strip().split('\t')
        if(flag == 0):
            ids.append(spl[0])
            length.append(spl[1])
        su = int(spl[2])+int(spl[3])
        co.append(str(su))
    count.append(co)
    co = []
    flag = 1
names = '\t'.join(fname)
print '\t'.join(['contig','length',names])
for i in range(0,len(ids)):
    tco = []
    for j in range(0,len(count)):
        tco.append(count[j][i])
    result =[ids[i],length[i],'\t'.join(tco)]
    print '\t'.join(result)
```

and running it with:

```bash
python get_count_table.py *.idxstats.txt > count.txt
```

And, finally, reducing the dimensions of the file by removing all contigs with less than 100 counts across all samples.

```bash
awk 'FNR == 1{print;next}{for(i=3;i<=NF;i++) if($i > 100){print;next}}' count.txt > count2.txt
```

