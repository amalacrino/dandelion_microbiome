### Run Kraken

```bash
DATADIR=$add_path_here #path to folder with fastq files after removing hsot plant reads
DBDIR=$add_path_here #path to the kraken database (see kraken2 manual)
OUTDIR=$add_path_here #path to kraken output directory
OUTDIR_bracken=$add_path_here #path to bracken output directory

cd $DATADIR

for file in *_host_removed_R1.fastq.gz
do
    SAMPLE=`basename $file _host_removed_R1.fastq.gz`
    cd $OUTDIR
    kraken2 --use-names --threads 64 \
      --db $DBDIR \
      --gzip-compressed \
      --paired $DATADIR/${SAMPLE}_host_removed_R1.fastq.gz $DATADIR/${SAMPLE}_host_removed_R2.fastq.gz \
      --output ${SAMPLE}.txt --report ${SAMPLE}.report.txt
done 

cd $OUTDIR

for file in *.report.txt
do
  SAMPLE=`basename $file .report.txt`
	bracken -d $DBDIR -i ${SAMPLE}.report.txt -o $OUTDIR_bracken/${SAMPLE}.txt -l S
done 
```

The content of `$OUTDIR_bracken` can be now downloaded using an SFTP client and processed locally using `R`. 

Before processing with `R`, place all the `.txt` files generated from bracken inside a folder called `data` and in the same directory place a file called metadata, which first column must report the sample ID that have so match the `.txt` filenames.

```R
library("phyloseq")

files <- dir(path = "data", pattern = "*.txt")
data <- sapply(files, function(x){
  raw.data <- read.table(paste0("data/", x), sep="\t", header = T)[,c("name", "new_est_reads")]
  colnames(raw.data)[2] = x
  return(raw.data)},
  simplify = FALSE,USE.NAMES = TRUE)
data <- Reduce(function(x,y) {merge(x,y,all=T, by="name")}, data)
data[is.na(data)]<-0
colnames(data) <- sub('\\.[^.]+$', '', colnames(data))

row.names(data) <- paste('Taxon', 1:nrow(data), sep = '_')
taxa <- as.data.frame(data[,c(1)])
row.names(taxa) <- row.names(data) 
colnames(taxa)[1] <- "Genus"
write.table(taxa, "Taxa.txt")

metadata <- read.table("metadata.txt", header = T, sep = "\t", row.names = 1)
metadata <- sample_data(metadata)

otutable <- as.matrix(data[,c(-1)])
otutable <- otu_table(otutable, taxa_are_rows=TRUE)

colnames(otutable) %in% rownames(metadata)

GM <- merge_phyloseq(otutable, metadata)

GM <- filter_taxa(GM, function (x) {sum(x > 0) > 1}, prune=TRUE)
GM <- subset_samples(GM, sample_type !="negctr_tapwater")
GM


GM.kracken <- GM
gdata::keep(GM.kracken, sure = TRUE)
save.image(file = 'ASV_table_kraken.rds')
```

This script generates the file `ASV_table_kraken.rds` that can be then used for further analyses. It also generates a file named `Taxa.txt` where it stores the correct taxonomy assignment along the dummy taxonomy.
