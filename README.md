# MEMEsuite

- 0.创建conda环境用于meme分析  
- 1.构建meme-chip所需的bed文件
- 2.meme-chip analysis   
- 3.Fimo analysis  

---

## 0.创建conda环境用于meme分析  

- Firstly, we create conda source to perform fimo analysis.  
- And then, we download the motif database from https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme.  
- We download the tf files from https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv.  

```bash
conda create -n meme
conda activate meme
conda install -c bioconda meme
conda install -c bioconda bedtools

# 生成一个filenames的文件，用来记录输出的文件名称（样本名称），例如：
ls bed500/*bed |cut -d "_" -f 2 |cut -d "/" -f 2 > filenames 
```

## 1.构建meme-chip所需的bed文件   
前面教程介绍过如何创建了，这里不再叙述  
- https://github.com/y741269430/featurecounts?tab=readme-ov-file#31-%E6%9E%84%E5%BB%BAmeme-chip%E6%89%80%E9%9C%80%E7%9A%84bed%E6%96%87%E4%BB%B6

## 2.meme-chip analysis   

将bed文件转为fasta文件  
```bash
vim meme_bed2fa.sh

#!/bin/bash
## BED to fa ##

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

# 设置参考基因组路径
ucsc_fa=/home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/GRCm39.genome.fa

# 读取文件名列表
cat filenames | while read i; do
    # 构造输入和输出文件路径
    input_file="$input_dir/${i}_equal_p.bed"
    output_file="$output_dir/${i}_mm10"

    # 运行 bedtools getfasta
    nohup bedtools getfasta -fi "$ucsc_fa" -bed "$input_file" -fo "$output_file" &
done

```
进行meme-chip分析  
```bash
vim memechip.sh

#!/bin/bash
## meme-chip ##

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

# 设置 motif 数据库路径
memedbs=/home/jjyang/downloads/Motif_database/merge_HM_JAS.meme

# 读取文件名列表
cat filenames | while read i; do
    # 构造输入和输出文件路径
    input_file="$input_dir/${i}_mm10"
    output_path="$output_dir/${i}_meme"

    # 运行 meme-chip
    nohup meme-chip -meme-p 6 -oc "$output_path" "$input_file" -db "$memedbs" &
done
```
    
## 3.Fimo analysis  

提前创建fimo文件夹
```bash
vim f1_fimo.sh

#!/bin/bash
## fimo ##

# 森检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

# 设置 motif 数据库路径
memedbs=/home/jjyang/downloads/Motif_database/merge_HM_JAS.meme

# 读取文件名列表
cat filenames | while read i; do
    # 构造输入和输出文件路径
    input_file="$input_dir/${i}_mm10"
    output_path="$output_dir/${i}"

    # 运行 fimo
    nohup fimo -oc "$output_path" "$memedbs" "$input_file" &
done
```

提取fimo文件夹中的 tsv 转换为 BED 进行peak注释（R）   
```bash
vim f2_tsv2bed.sh

#!/bin/bash
## tsv2bed ##

path=./bed500

cat filenames | while read i; 
do  
cat ./fimo/${i}/fimo.tsv |awk 'NR ==1 {next} {print $2"\t"$3"\t"$4"\t"$1"\t"$7}' |awk '/chr/ {print $1"\t"strtonum($2)"\t"strtonum($3)"\t"$4"\t"$5}' > $path/${i}_fimo.bed &
done
```

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



