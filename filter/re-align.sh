#! /bin/bash

echo "Initializing..."
export BASE_PATH=/mnt/raid6/bacphagenetwork/data/tmp
export INDEXING_PATH=/mnt/raid6/bacphagenetwork/data/00_bwa_index/chm13v2/chm13v2.0_noY.fa
export GENE_DATA_1=${BASE_PATH}/BJ001_1.fastq.gz
export GENE_DATA_2=${BASE_PATH}/BJ001_2.fastq.gz

echo "Re-aligning..."
bwa mem -t 4 -p $INDEXING_PATH $GENE_DATA_1 $GENE_DATA_2 -M > ${BASE_PATH}/BJ001.sam
echo "Re-alignment done."

echo "===================="

echo "Filtering..."
samtools view -b -F 4 ${BASE_PATH}/BJ001.sam > ${BASE_PATH}/BJ001.bam
echo "Filtering done."

echo "===================="

echo "Sorting and indexing..."
samtools sort ${BASE_PATH}/BJ001.bam -o ${BASE_PATH}/BJ001.sorted.bam
samtools index ${BASE_PATH}/BJ001.sorted.bam

echo "===================="

echo "Format converting..."
samtools view -h ${BASE_PATH}/BJ001.sorted.bam > ${BASE_PATH}/BJ001.sorted.sam
echo "Format converting done."

echo "===================="

echo "All done."
