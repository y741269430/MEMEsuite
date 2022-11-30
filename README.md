# MEMEsuite

Firstly, we create conda source to perform fimo analysis.  
And then, we download the motif database from https://meme-suite.org/meme/doc/download.html.  

    conda create -n meme
    conda activate meme
    conda install -c bioconda meme
    conda install -c bioconda bedtools

    ls pm_saf/*bed |cut -d "_" -f 2 |cut -d "/" -f 2 > filenames

## BED files  

Firstly, we should make the BED files for downstream analysis. The BED files were make in R.  
For example:  
We annotated the BED files by ChIPseeker and ChIPpeakAnno.  

    library(ChIPseeker)
    library(ChIPpeakAnno)

    peak <- lapply(list.files('ATAC-nt-rawdata/peak/', "*.bed"), 
                   function(x){return(readPeakFile(file.path('ATAC-nt-rawdata/peak/', x)))})
    names(peak) <- c('e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5')
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 0), 
                           annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
    peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})
    
And then we divided the chromatin into promoter regions and gene body regions.  
    
    region_bed <- lapply(peakAnno_df, function(x){
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
    pm_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
        x <- x[grep("Promoter", ignore.case = F, x$annotation), ]
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
    gb_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
        x <- x[c(grep("5' UTR", ignore.case = F, x$annotation),
                 grep("Intron", ignore.case = F, x$annotation),
                 grep("Exon", ignore.case = F, x$annotation),
                 grep("Downstream", ignore.case = F, x$annotation),
                 grep("3' UTR", ignore.case = F, x$annotation)), ]
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
    dis_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL), ]
        x <- x[grep("Distal Intergenic", ignore.case = F, x$annotation),]
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
  
Finally, we saved the regions files to BED.  

    for (i in 1:length(region_bed)) {
        write.table(region_bed[i],
                    paste('ATAC-nt-rawdata/saf/', names(region_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }
    for (i in 1:length(pm_bed)) {
        write.table(pm_bed[i],
                    paste('ATAC-nt-rawdata/pm_saf/', names(pm_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }
    for (i in 1:length(gb_bed)) {
        write.table(gb_bed[i],
                    paste('ATAC-nt-rawdata/gb_saf/', names(gb_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }
    for (i in 1:length(dis_bed)) {
        write.table(dis_bed[i],
                    paste('ATAC-nt-rawdata/dis_saf/', names(dis_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }

The BED files can be used to convert to saf files for featurecount, it also can be used to fimo analysis.  

----

## Make BED files to equal length BED files  

    vim f1_bed2equal.sh

    #!/bin/bash
    ## make BED2equal.config ##
    ## BED to equal length BED ##

    cat filenames | while read i; 
    do
    nohup awk -v FS="\t" -v OFS="\t" '{midpos=$2+$5;print $1,midpos-250,midpos+250;}' ./pm_saf/${i}_allpeak.bed > ./pm_saf/${i}_equal_p.bed &
    done

## Convert BED files to fa files   

    vim f2_bed2fa.sh

    #!/bin/bash
    ## BED to fa ##

    ucsc_fa=/home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/GRCm38.primary_assembly.genome.fa

    cat filenames | while read i; 
    do   
    bedtools getfasta -fi $ucsc_fa -bed ./pm_saf/${i}_equal_p.bed -fo ./pm_saf/${i}_mm10 &
    done

## Fimo analysis     

    vim f3_fimo.sh

    #!/bin/bash
    ## fimo ##

    memedbs=/home/yangjiajun/downloads/Motif_database/merge_HM_JAS.meme

    cat filenames | while read i; 
    do
    nohup fimo -oc ./pm_saf/${i} $memedbs ./pm_saf/${i}_mm10 &
    done

## Extract the tsv files and convert to BED files

    vim f4_tsv2bed.sh

    #!/bin/bash
    ## tsv2bed ##

    cat filenames | while read i; 
    do  
    cat ./pm_saf/${i}/fimo.tsv |awk 'NR ==1 {next} {print $2"\t"$1"\t"$7}' |awk '{gsub(/:|-/, "\t", $1); print $0}' |awk '!a[$0]++{print}' |awk '/chr/ {print $1"\t"strtonum($2)"\t"strtonum($3)"\t"$4"\t"$5}' > ./pm_saf/${i}_fimo.bed &
    done
