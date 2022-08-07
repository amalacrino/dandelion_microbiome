As first step, we need to assemble the row reads into contigs (basically we merge short reads into a longer read). For this step we use MegaHit. From this experiment we obtained a massive amount of data, and MegaHit requires too much memory to process it all together. So, here, we are going to perform the contig assebly separately for each sample and later in the pipeline we will merge this output. Using the raw reads filtered from the host reads as input, we can run this code:

```bash
DATADIR=$_add_path_here
OUTDIR=$_add_path_here

cd $DATADIR

for file in *_host_removed_R1.fastq.gz
do
    SAMPLE=`basename $file _host_removed_R1.fastq.gz`
		megahit -t 32 -1 ${SAMPLE}_host_removed_R1.fastq.gz \
			-2 ${SAMPLE}_host_removed_R2.fastq.gz \
			-o $OUTDIR/${SAMPLE}
done 
```

MegaHit creates a folder for each of our inputs, with a file named `final.contigs.fa` were it stores the information we will need for the nexrt steps. So, running the code below we will create a new folder, take each `final.contigs.fa` file, rename it according to the folder name (which is our sample), and move it to the new folder.

```bash
MegahitPATH=$_add_path_here
MegahitFILES=$_add_path_here
find "$MegahitPATH" -type f -iname 'final.contigs.fa' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" "$MegahitFILES/${path##*/}.fa" ' sh_cp {} \;

cd $MegahitFILES
cat *.fa > allseqs.fa
```

Here is a tip. MegaHit gives the same names to the contigs obtained from different files (something like contig1. contig2, ...) and some softwares down in the pipeline do not like to have contigs with the same name within the same `fasta` file. So, here we are going to rename all the contigs in the file we just obtained. First, we need to know how many contigs our file has:

```bash
cd 4c_mh_all
sed '/^>/d' allseqs.fa | wc -l
```
In our case the result was `29853253`. So we are going to use this number to run the line below.

```bash
for i in {1..29853253}; do echo A_$i; done | paste - <(sed '/^>/d' allseqs.fa) | sed -e 's/^/>/' -e 's/\t/\n/' > new_file.fa
```


Remove duplicate reads

```bash
#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem=16G
#SBATCH --time=00:10:00
#SBATCH --job-name=removedups
#SBATCH --mail-type=ALL --mail-user=antonino.malacrino@gmail.com
#SBATCH --account=pas1910

PRJDIR=/users/PAS1609/antoninomalacrino/bennett/selenium_metagenomics
DATADIR=$PRJDIR/4_megahit2

module load python/3.7-2019.10
source activate $HOME/conda_envs/seqkit

cd $DATADIR

seqkit rmdup -s < new_file.fa > input.faa

rm allseqs.fa new_file.fa
```

Run Prokka:

```bash
#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=32 --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=runPROKKA
#SBATCH --mail-type=ALL --mail-user=antonino.malacrino@gmail.com
#SBATCH --account=pas1910

PRJDIR=/users/PAS1609/antoninomalacrino/bennett/selenium_metagenomics
DATADIR=$PRJDIR/4_megahit2
OUTDIR=$PRJDIR/5_prokka

cp $DATADIR/input.faa $OUTDIR/input.faa 

cd $OUTDIR

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < input.faa

source activate $HOME/conda_envs/prokka

find -name "*.fa" | parallel -j 32 prokka --outdir $OUTDIR --norrna --notrna --metagenome --force --cpus 1 --prefix {}\.prokka {}
```

Collate tsv

```bash
mkdir 6_prokka2
cd 5_prokka
INDIR=/users/PAS1609/antoninomalacrino/bennett/selenium_metagenomics/6_prokka2
cp *.tsv $INDIR
cd ../6_prokka2
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.tsv > prokka.txt
rm *.tsv
```

Run bowtie

```bash
mkdir 7_bowtie
mkdir 7_bowtie/input 7_bowtie/output 7_bowtie/final_output
```

```bash
#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=16 --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=runBowtie
#SBATCH --mail-type=ALL --mail-user=antonino.malacrino@gmail.com
#SBATCH --account=pas1910

module load python/3.7-2019.10
source activate $HOME/conda_envs/bowtie2

PRJDIR=/users/PAS1609/antoninomalacrino/bennett/selenium_metagenomics

DATADIR=$PRJDIR/5_prokka
RAWREADS=$PRJDIR/2_trimgalore
INDIR=$PRJDIR/7_bowtie/input
OUTDIR=$PRJDIR/7_bowtie/output
OUTDIR2=$PRJDIR/7_bowtie/final_output

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

conda deactivate

source activate $HOME/conda_envs/samtools

for x in *.sam
do
  f=`basename $x .sam`;
  echo $f;
  samtools view -@ 16 -bS ${f}.sam > ${f}.bam;
  samtools sort -@ 16 ${f}.bam -o ${f}.sorted;
  samtools index -@ 16 ${f}.sorted;
  samtools idxstats ${f}.sorted > ${f}.idxstats.txt;
done

```

Save this file as `get_count_table.py`

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

Run

```bash
#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --mem=16G
#SBATCH --time=06:00:00
#SBATCH --job-name=getCounts
#SBATCH --mail-type=ALL --mail-user=antonino.malacrino@gmail.com
#SBATCH --account=pas1910

PRJDIR=/users/PAS1609/antoninomalacrino/bennett/selenium_metagenomics
DATADIR=$PRJDIR/7_bowtie/final_output

cd $DATADIR

module load python/2.7-conda5.2

python get_count_table.py *.idxstats.txt > count.txt

awk 'FNR == 1{print;next}{for(i=3;i<=NF;i++) if($i > 100){print;next}}' count.txt > count2.txt
```






