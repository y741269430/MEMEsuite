# MEMEsuite

We download the docker image from https://hub.docker.com/r/memesuite/memesuite  

    docker pull memesuite/memesuite:latest

conda create -n meme
conda activate meme
conda install -c bioconda meme

ls pm_saf/*bed |cut -d "_" -f 2 |cut -d "/" -f 2 > filenames

vim f1_bed2equal.sh

#!/bin/bash
## make bed2equal.config ##
## region bed to equal length bed ##

ls ~/ATAC-nt-rawdata/saf/*_allpeak.bed > 1
ls ~/ATAC-nt-rawdata/saf/*_allpeak.bed |cut -d "/" -f 6 |cut -d "." -f 1 > 0
paste 0 1 > bed2equal.config

cat bed2equal.config | while read id;
do echo $id
arr=($id)
id1=${arr[1]}
sample=${arr[0]}

nohup awk -v FS="\t" -v OFS="\t" '{midpos=$2+$5;print $1,midpos-250,midpos+250;}' $id1 > ~/ATAC-nt-rawdata/memefimo/${sample}/${sample}_equal_p.bed &
done

rm 0 1 bed2equal.config

