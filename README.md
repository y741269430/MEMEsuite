# MEMEsuite

- 0.创建conda环境用于meme分析  
- 1.构建meme-chip所需的bed文件
- 2.peak取重叠区域做meme分析
- 3.bed转fasta（巨坑）
- 4.meme-chip analysis   
- 5.Fimo analysis  

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
ls bed/*bed |cut -d "_" -f 2 |cut -d "/" -f 2 > filenames 
```

## 1.构建meme-chip所需的bed文件   
前面教程介绍过如何创建了，这里不再叙述  
- [featurecounts.md](https://github.com/y741269430/featurecounts?tab=readme-ov-file#31-%E6%9E%84%E5%BB%BAmeme-chip%E6%89%80%E9%9C%80%E7%9A%84bed%E6%96%87%E4%BB%B6)

## 2.peak取重叠区域做meme分析
前面教程讲了如何构建重叠区域  
- [计算重复样本的 peak 之间的重叠的坐标位置.md](https://github.com/y741269430/ATAC-seq/blob/main/%E8%AE%A1%E7%AE%97%E9%87%8D%E5%A4%8D%E6%A0%B7%E6%9C%AC%E7%9A%84%20peak%20%E4%B9%8B%E9%97%B4%E7%9A%84%E9%87%8D%E5%8F%A0%E7%9A%84%E5%9D%90%E6%A0%87%E4%BD%8D%E7%BD%AE.md)

```bash
bedtools intersect -a CON_1.bed -b CON_2.bed > intersect_CON.bed
bedtools intersect -a treatment_1.bed -b treatment_2.bed > intersect_treatment.bed
```

> 这里会有错误：在处理基因组peak数据时，如果对重叠区域进行筛选后，每个peak的长度变得不一致，这将导致MEME工具无法正确识别并输出motif，最终得到空的结果。为了解决这一问题，我们可以采取以下两种策略之一：    
> 1.独立样本输入：放弃对所有样本取交集以获得共同的peak区域的做法，转而将每个样本的原始peak数据分别作为MEME的输入。这种方法保留了每个样本的独特信息，避免了因长度不一致而导致的问题。    
> 2.统一peak长度：在构建重叠区域时，确保所有peak具有相同的长度。可以通过将较短的peak扩展至一个固定的长度来实现，例如，将peak的末端位置调整为其起始端加300碱基对（bp）。这样做的目的是使所有输入到MEME的peak拥有相同长度，从而保证分析的有效性。    

```bash
bedtools intersect -a CON_1.bed -b CON_2.bed | awk -F'\t' '{print $1, $2, $2+300, $4}' OFS='\t' > intersect_CON.bed
bedtools intersect -a treatment_1.bed -b treatment_2.bed | awk -F'\t' '{print $1, $2, $2+300, $4}' OFS='\t' > intersect_treatment.bed
```

## 3.bed转fasta（巨坑）   
参考   
- [BED Genomic Loci Format](https://meme-suite.org/meme/doc/bed-format.html)  
- [BED2FASTA](https://meme-suite.org/meme/doc/bed2fasta.html?man_type=web)  

使用以下两个命令处理基因组数据时，会得到不同的输出结果，meme需要输入的是第一个命令产生的结果。而bedtools getfasta命令产生的结果文件会有error，输出空值。因此推荐使用以下meme_bed2fa.sh脚本    
```bash
bed2fasta -name intersect_CON.bed /home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/GRCm39.genome.fa > test
bedtools getfasta -name -bed intersect_CON.bed -fi /home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/GRCm39.genome.fa -fo test2
```
<img src="https://github.com/y741269430/MEMEsuite/blob/main/img/bed2fa%E7%BB%93%E6%9E%9C%E5%AF%B9%E6%AF%94.png" width="600" />

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
    input_file="$input_dir/intersect_${i}.bed"
    output_file="$output_dir/${i}_mm10"

    # 运行 bedtools getfasta
    nohup bed2fasta -name "$input_file" "$ucsc_fa" > "$output_file" &
done
```
执行该脚本  第一个是input 第二个是output  
```bash
conda activate meme    
bash meme_bed2fa.sh peak200/ peak200/
```

## 4.meme-chip analysis   
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
执行该脚本  第一个是input 第二个是output  
```bash
bash memechip.sh peak200/ meme_chip_result/
```

## meme输出结果的error  
<img src="https://github.com/y741269430/MEMEsuite/blob/main/img/meme_error.png" width="900" />   

> 这个错误大概意思是 Open MPI 在尝试加载 librdmacm.so.1 这个共享对象文件时失败了，原因是找不到这个文件，并且在使用 mpirun 运行 meme 程序时出现了问题，其中一个或多个进程以非零状态退出，导致整个任务被终止。此外，错误信息还提到了 CMA（Contiguous Memory Allocator）权限被拒绝的问题。   
初步解决方案      

```bash
sudo apt-get install libibverbs1

echo $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
source ~/.bashrc
```
> 当我重新运行memechip.sh脚本的时候，结果文件又有error了    
<img src="https://github.com/y741269430/MEMEsuite/blob/main/img/meme_error2.png" width="600" />

> 大概意思是有重复的基因名，它需要唯一的基因。。没办法只好把名字改成SYMBOL_峰值位置（这样总不会重复了吧）    
<img src="https://github.com/y741269430/MEMEsuite/blob/main/img/centrimo_error.png" width="600" />

> 至于centrimo 这个问题，当时计划取的是peak峰值左右100bp的区间，作为motif预测的区域，结果在bed转fasta的时候，有些碱基是被去掉的，导致大部分长度没到200bp，它默认给的参数是`-seqlen 200 `。但是这个命令是有问题的，当我取400bp放进去的时候，它就会报错说你不到400bp，它默认给的参数变为了`-seqlen 400 `。解决方案就是，单独跑这个centrimo 程序

## 5.Fimo analysis  

提前创建fimo文件夹
```bash
vim f1_fimo.sh

#!/bin/bash
## fimo ##

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
    output_path="$output_dir/${i}"

    # 运行 fimo
    nohup fimo -oc "$output_path" "$memedbs" "$input_file" &
done
```
执行该脚本  第一个是input 第二个是output  
```bash
bash f1_fimo.sh peak200/ fimo_results/
```
提取fimo文件夹中的 tsv 转换为 BED 进行peak注释（R）   
```bash
vim f2_tsv2bed.sh

#!/bin/bash
## tsv2bed ##

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 获取输入和输出目录路径
input_dir=$1
output_dir=$2

cat filenames | while read i; 
do  
input_file="$input_dir/${i}/fimo.tsv"
output_path="$output_dir/${i}_fimo.bed"

cat $input_file |awk 'NR ==1 {next} {print $2"\t"$3"\t"$4"\t"$1"\t"$7}' \
|awk '/chr/ {print $1"\t"strtonum($2)"\t"strtonum($3)"\t"$4"\t"$5}' > $output_path &
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



