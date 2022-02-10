#!bin/bash
rawData_path=/data2/project-NPC/RNA-seq/rawData
cleanData_path=/data2/project-NPC/RNA-seq/cleanData

cat id.txt | while read line 
do
  java -jar /home/Soft_ware/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
    -phred33 \
    -threads 36 \
    ${rawData_path}/${line}_1.fastq.gz \
    ${rawData_path}/${line}_2.fastq.gz \
    ${cleanData_path}/${line}_1.paired.fastq.gz \
    ${cleanData_path}/${line}_1.unpaired.fastq.gz \
    ${cleanData_path}/${line}_2.paired.fastq.gz \
    ${cleanData_path}/${line}_2.unpaired.fastq.gz \
    LEADING:5  TRAILING:5  SLIDINGWINDOW:4:20  MINLEN:50 HEADCROP:10 
done
