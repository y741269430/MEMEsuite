# MEMEsuite

- 0.Build source
- 1.BED files (R)  
- 2.Make BED files to equal length BED files 
- 3.Convert BED files to fa files  
- 4.Fimo analysis   
- 5.Extract the tsv files and convert to BED files
- 6.Read BED files and construct the connection files (R)    

----

## 0.Build source  

Firstly, we create conda source to perform fimo analysis.  
And then, we download the motif database from https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme.  
We download the tf files from https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv.  

    conda create -n meme
    conda activate meme
    conda install -c bioconda meme
    conda install -c bioconda bedtools

    ls pm_saf/*bed |cut -d "_" -f 2 |cut -d "/" -f 2 > filenames

## 1.BED files (R)  

Firstly, we should make the BED files for downstream analysis. The BED files were make in R.  
For example:  
We annotated the BED files by ChIPseeker and ChIPpeakAnno.  

    library(ChIPseeker)
    library(ChIPpeakAnno)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)

    peak <- lapply(list.files('ATAC-nt-rawdata/peak/', "*.bed"), 
                   function(x){return(readPeakFile(file.path('ATAC-nt-rawdata/peak/', x)))})
    names(peak) <- c('e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5')
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 0), 
                           annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
    peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})  
    
And then we divided the chromatin into promoter regions and gene body regions.  

![peak_anno.png](https://github.com/y741269430/MEMEsuite/blob/main/peak_anno.png)  

    region_bed <- lapply(peakAnno_df, function(x){
        # colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        # x <- x[, c(1,2,3,23,7,8)]
        x <- x[, c(1,2,3,19,7,5)]
        return(x)
      })
    pm_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
        x <- x[grep("Promoter", ignore.case = F, x$annotation), ]
        # colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        # x <- x[, c(1,2,3,23,7,8)]
        x <- x[, c(1,2,3,19,7,5)]
        return(x)
      })
    gb_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
        x <- x[c(grep("5' UTR", ignore.case = F, x$annotation),
                 grep("Intron", ignore.case = F, x$annotation),
                 grep("Exon", ignore.case = F, x$annotation),
                 grep("Downstream", ignore.case = F, x$annotation),
                 grep("3' UTR", ignore.case = F, x$annotation)), ]
        # colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        # x <- x[, c(1,2,3,23,7,8)]
        x <- x[, c(1,2,3,19,7,5)]
        return(x)
      })
    dis_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL), ]
        x <- x[grep("Distal Intergenic", ignore.case = F, x$annotation),]
        # colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        # x <- x[, c(1,2,3,23,7,8)]
        x <- x[, c(1,2,3,19,7,5)]
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

    vim bed2saf.sh

    #!/bin/bash
    ## make saf (bedtools) for featurecount ##

    path=./pm_saf
    
    cat filenames | while read i; 
    do
    nohup bedtools sort -i $path/${i}_allpeak.bed > $path/${i}.rm.bed && bedtools merge -c 4,6 -o first -i $path/${i}.rm.bed |awk 'BEGIN{print "GeneID" "\t"  "Chr" "\t" "Start" "\t" "End" "\t" "Strand"}{print $4"\t"$1"\t"strtonum($2)"\t"strtonum($3)"\t"$5}' > $path/${i}.saf && rm $path/${i}.rm.bed -rf &
    done

----

## 2.Make BED files to equal length BED files  

    vim f1_bed2equal.sh

    #!/bin/bash
    ## make BED2equal.config ##
    ## BED to equal length BED ##
    
    path=./pm_saf
    
    cat filenames | while read i; 
    do
    nohup awk -v FS="\t" -v OFS="\t" '{midpos=$2+$5;print $1,midpos-250,midpos+250;}' $path/${i}_allpeak.bed > $path/${i}_equal_p.bed &
    done

## 3.Convert BED files to fa files   

    vim f2_bed2fa.sh

    #!/bin/bash
    ## BED to fa ##
    
    path=./pm_saf
    ucsc_fa=/home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/GRCm38.primary_assembly.genome.fa

    cat filenames | while read i; 
    do   
    bedtools getfasta -fi $ucsc_fa -bed $path/${i}_equal_p.bed -fo $path/${i}_mm10 &
    done

## 4.Fimo analysis     

    vim f3_fimo.sh

    #!/bin/bash
    ## fimo ##
    
    path=./pm_saf
    memedbs=/home/yangjiajun/downloads/Motif_database/merge_HM_JAS.meme

    cat filenames | while read i; 
    do
    nohup fimo -oc $path/${i} $memedbs $path/${i}_mm10 &
    done

## 5.Extract the tsv files and convert to BED files  

    vim f4_tsv2bed.sh

    #!/bin/bash
    ## tsv2bed ##

    path=./pm_saf
    
    cat filenames | while read i; 
    do  
    cat $path/${i}/fimo.tsv |awk 'NR ==1 {next} {print $2"\t"$1"\t"$7}' |awk '{gsub(/:|-/, "\t", $1); print $0}' |awk '!a[$0]++{print}' |awk '/chr/ {print $1"\t"strtonum($2)"\t"strtonum($3)"\t"$4"\t"$5}' > $path/${i}_fimo.bed &
    done

## 6.Read BED files and construct the connection files (R)    

Firstly, we read and annotated the fimo.bed files generated above.  

    library(ChIPseeker)
    library(ChIPpeakAnno)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)

    path <- 'pm_saf/'
    peak <- lapply(list.files(paste0('ATAC-nt-rawdata/fimo/', path), "*_fimo.bed"), function(x){
      return(readPeakFile(file.path(paste0('ATAC-nt-rawdata/fimo/', path), x)))
    })
    names(peak) <- c('e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5')

    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 0), 
                           annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
    peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})

    save(peakAnno_df, file = paste0('ATAC-nt-rawdata/fimo/', path, 'nt_fimo_Anno.RData'))

Download and read the tf files from https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv  

    tf <- read.csv("downloads/Motif_database/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv", sep = '\t')[,c(1,2,17,4,8,14)]
    colnames(tf)
    [1] "Model"                "Transcription.factor" "EntrezGene"           "Quality"              "Data.source"          "TF.family"

Merge peakAnno_df and tf  

![fimo_peak_anno.png](https://github.com/y741269430/MEMEsuite/blob/main/fimo_peak_anno.png)  

    fimo_link <- lapply(peakAnno_df, function(x){
      x <- x[,c(1,2,3,6,7,8,12,14,16:19)]
      colnames(x) <- c("seqnames", "start", "end", "Model", "pvalue", "annotation", "geneLength", "ENTREZID", "distanceToTSS", "ENSEMBL", "targetGene", "GENENAME")
      x <- merge(tf, x, "Model")
      return(x)
    })

We selected transcription factors validated by ChiP Seq and codding genes.  

    fimo_linkcod <- lapply(fimo_link, function(x){
      x <- x[-grep("cDNA", ignore.case = T, x$GENENAME),]
      x <- x[-grep("pseudogene", ignore.case = T, x$GENENAME),]
      x <- x[-grep("small nucleolar RNA", ignore.case = T, x$GENENAME),]
      x <- x[-grep("predicted gene", ignore.case = T, x$GENENAME),]
      x <- x[grep("ChIP-Seq", ignore.case = T, x$Data.source),]
      return(x)
    })

We connect the transcription factors and their target genes.  

    link_genes <- lapply(fimo_linkcod, function(x){x <- paste(x$Transcription.factor, x$targetGene, sep = "~")})

    co <- Reduce(intersect, list(link_genes[[1]],
                                 link_genes[[2]],
                                 link_genes[[3]],
                                 link_genes[[4]],
                                 link_genes[[5]]))
    co <- unique(co)
    co <- stringr::str_split_fixed(co, "~", n = 2)

    fimo_link_merge <- data.frame(regulatoryGene = co[,1], targetGene = co[,2])
    fimo_link_each <- lapply(fimo_linkcod, function(x){x <- data.frame(regulatoryGene = x[,2], targetGene = x[,16])})
    
    write.csv(fimo_link_merge, paste0('ATAC-nt-rawdata/fimo/', path, 'fimo_link_merge.csv'), row.names = F)
    save(fimo_link, fimo_linkcod, fimo_link_merge, fimo_link_each, file = paste0('ATAC-nt-rawdata/fimo/', path, 'nt_fimo_link.Rdata'))




