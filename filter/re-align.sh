#! /bin/bash

echo "Initializing..."
export BASE_PATH=/mnt/raid6/bacphagenetwork/data/tmp
export WORKING_PATH=/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/filter
export INDEXING_PATH=/mnt/raid6/bacphagenetwork/data/00_bwa_index/chm13v2/chm13v2.0_noY.fa
export GENE_DATA_1=${BASE_PATH}/BJ001_1.fastq.gz
export GENE_DATA_2=${BASE_PATH}/BJ001_2.fastq.gz

echo "Re-aligning..."
bwa mem -t 16 $INDEXING_PATH $GENE_DATA_1 $GENE_DATA_2 -M > ${BASE_PATH}/BJ001.sam 2> ${WORKING_PATH}/bwa_mem.log
echo "Re-alignment done."

echo "===================="

echo "Filtering..."
samtools view -b -F 4 ${BASE_PATH}/BJ001.sam > ${BASE_PATH}/BJ001.bam 2> ${WORKING_PATH}/samtools_view.log
echo "Filtering done."

echo "===================="

echo "Sorting and indexing..."
samtools sort ${BASE_PATH}/BJ001.bam -o ${BASE_PATH}/BJ001.sorted.bam 2> ${WORKING_PATH}/samtools_sort.log
samtools index ${BASE_PATH}/BJ001.sorted.bam

echo "===================="

echo "Format converting..."
samtools view -h ${BASE_PATH}/BJ001.sorted.bam > ${BASE_PATH}/BJ001.sorted.sam 2> ${WORKING_PATH}/samtools_view2.log
echo "Format converting done."

echo "===================="

echo "All done."
