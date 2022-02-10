#!bin/bash
clean_path=/data2/project-NPC/RNA-seq/cleanData
genome_path=/data/DataBase/genomes/humman/ENSEMBLE
mapping_path=/data2/project-NPC/RNA-seq/mapping
# build a hisat2 index
hisat2_extract_splice_sites.py ${genome_path}/Homo_sapiens.GRCh38.100.gtf >Homo_sapiens.GRCh38.ss
hisat2_extract_exons.py ${genome_path}/Homo_sapiens.GRCh38.100.gtf >Homo_sapiens.GRCh38.exon
hisat2-build --ss Homo_sapiens.GRCh38.ss --exon Homo_sapiens.GRCh38.exon ${genome_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38_tran

cat id.txt | while read line

do
	hisat2 -p 32 --dta -x ${genome_path}/Homo_sapiens.GRCh38_tran -1 ${clean_path}/${line}_1.paired.fastq.gz -2 ${clean_path}/${line}_2.paired.fastq.gz -S ${mapping_path}/${line}.sam
	samtools view -S -@ 16 ${mapping_path}/${line}.sam -b > ${mapping_path}/${line}.bam 
  samtools sort -@ 16 ${mapping_path}/${line}.bam -o ${mapping_path}/${line}_sorted.bam
  rm ${mapping_path}/${line}.sam ${mapping_path}/${line}.bam
done

# Gene expression was quantified using featureCounts
featureCounts -T 32 -p -t exon -g gene_id -a ${genome_path}/Homo_sapiens.GRCh38.100.gtf -o ${mapping_path}/GSE102349_featurecounts.txt ${mapping_path}/*_sorted.bam
