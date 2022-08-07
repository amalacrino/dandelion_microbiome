We used VSEARCH to process raw data using the code below:

```bash
PRJ_DIR=$add_path_here
RAW_DATA=$add_path_here #this path is different for 16S and ITS data
REFERENCE=$add_path_here #this path is reference.fa for 16S or unite.fa for ITS
THREADS=16
CLUSTERID=0.97
MAXEE=1.0

cd $REFERENCE
wget https://www.drive5.com/sintax/silva_16s_v123.fa.gz
gunzip silva_16s_v123.fa.gz
mv silva_16s_v123.fa reference.fa

wget https://files.plutof.ut.ee/public/orig/82/CB/82CB44BBAAA7D3AEAC297B5689BDA2963E8D0666E01FE0B54096147AFAF85263.gz
gunzip 82CB44BBAAA7D3AEAC297B5689BDA2963E8D0666E01FE0B54096147AFAF85263.gz
mv 82CB44BBAAA7D3AEAC297B5689BDA2963E8D0666E01FE0B54096147AFAF85263 unite.fa

# Enter subdirectory

cd $PRJ_DIR
mkdir $PRJ_DIR/VSEARCH/
mkdir $PRJ_DIR/VSEARCH/1_QC
mkdir $PRJ_DIR/VSEARCH/1_DereplicateSS
mkdir $PRJ_DIR/VSEARCH/1_DereplicateAS
mkdir $PRJ_DIR/VSEARCH/2_Cluster

# Process samples

for f in $RAW_DATA/*.fastq; do

    s=`basename $f .fastq`;

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================
    
    echo
    echo Quality filtering
    echo

    $VSEARCH --fastq_filter $RAW_DATA/$s.fastq \
        --threads $THREADS \
        --fastq_maxee $MAXEE \
        --fastq_minlen 100 \
        --fastq_maxlen 600 \
        --fastq_maxns 50 \
        --fastaout $PRJ_DIR/VSEARCH/1_QC/$s.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample.n
    echo

    $VSEARCH --derep_fulllength $PRJ_DIR/VSEARCH/1_QC/$s.filtered.fasta \
        --threads $THREADS \
        --strand plus \
        --output $PRJ_DIR/VSEARCH/1_DereplicateSS/$s.derep.fasta \
        --sizeout \
        --relabel $s. \
        --fasta_width 0

done

echo
echo ====================================
echo Processing all samples together
echo ====================================
echo
echo Merge all samples

cat $PRJ_DIR/VSEARCH/1_DereplicateSS/*.derep.fasta > $PRJ_DIR/VSEARCH/1_DereplicateAS/all.fasta

echo
echo Sum of unique sequences in each sample: $(cat $PRJ_DIR/VSEARCH/1_DereplicateAS/all.fasta | grep -c "^>")
echo
echo Dereplicate across samples
echo

$VSEARCH --derep_fulllength $PRJ_DIR/VSEARCH/1_DereplicateAS/all.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc $PRJ_DIR/VSEARCH/1_DereplicateAS/all.derep.uc \
    --output $PRJ_DIR/VSEARCH/1_DereplicateAS/derep.fasta

echo
echo Unique non-singleton sequences: $(grep -c "^>" $PRJ_DIR/VSEARCH/1_DereplicateAS/derep.fasta)
echo
echo Cluster sequences using VSEARCH
echo

$VSEARCH --cluster_size $PRJ_DIR/VSEARCH/1_DereplicateAS/derep.fasta \
    --threads $THREADS \
    --id $CLUSTERID \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --centroids $PRJ_DIR/VSEARCH/2_Cluster/centroids.fasta

echo
echo Clusters: $(grep -c "^>" $PRJ_DIR/VSEARCH/2_Cluster/centroids.fasta)
echo
echo Sort and remove singletons
echo

$VSEARCH --sortbysize $PRJ_DIR/VSEARCH/2_Cluster/centroids.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --minsize 2 \
    --output $PRJ_DIR/VSEARCH/2_Cluster/sorted.fasta

echo
echo Non-singleton clusters: $(grep -c "^>" $PRJ_DIR/VSEARCH/2_Cluster/sorted.fasta)
echo 
echo De novo chimera detection
echo

$VSEARCH --uchime_denovo $PRJ_DIR/VSEARCH/2_Cluster/sorted.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --nonchimeras $PRJ_DIR/VSEARCH/2_Cluster/nonchimeras.fasta \

echo
echo Unique sequences after de novo chimera detection: $(grep -c "^>" $PRJ_DIR/VSEARCH/2_Cluster/nonchimeras.fasta)

echo
echo Relabel OTUs
echo

$VSEARCH --fastx_filter $PRJ_DIR/VSEARCH/2_Cluster/nonchimeras.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --relabel OTU_ \
    --fastaout $PRJ_DIR/VSEARCH/2_Cluster/otus.fasta

echo
echo Number of OTUs: $(grep -c "^>" $PRJ_DIR/VSEARCH/2_Cluster/otus.fasta)
echo
echo Map sequences to OTUs by searching
echo

$VSEARCH --usearch_global $PRJ_DIR/VSEARCH/1_DereplicateAS/all.fasta \
    --threads $THREADS \
    --db $PRJ_DIR/VSEARCH/2_Cluster/otus.fasta \
    --id $CLUSTERID \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --dbmask none \
    --otutabout $PRJ_DIR/VSEARCH/2_Cluster/otutab.txt

echo
echo Sort OTU table numerically
echo

sort -k1.5n $PRJ_DIR/VSEARCH/2_Cluster/otutab.txt > $PRJ_DIR/VSEARCH/2_Cluster/otutab.sorted.txt

echo
echo Done

mkdir $PRJ_DIR/VSEARCH/3_Taxonomy

vsearch -sintax $PRJ_DIR/VSEARCH/2_Cluster/otus.fasta -db $REFERENCE -tabbedout $PRJ_DIR/VSEARCH/3_Taxonomy/ASV_tax_raw.txt -strand both -sintax_cutoff 0.5 -threads $THREADS
```
