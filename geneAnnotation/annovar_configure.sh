#!/bin/bash

# 设置环境变量
ANNOVAR_DIR="/mnt/raid6/bacphagenetwork/tools//annovar"
DB_DIR="${ANNOVAR_DIR}/humandb"  # 数据库目录
BUILD="hg38"  # 基因组版本

# 创建数据库目录
mkdir -p "${DB_DIR}"

echo ">>> 下载RefGene数据库..."
perl "${ANNOVAR_DIR}/annotate_variation.pl" -buildver "${BUILD}" -downdb -webfrom annovar refGene "${DB_DIR}"

echo ">>> 下载ClinVar数据库..."
perl "${ANNOVAR_DIR}/annotate_variation.pl" -buildver "${BUILD}" -downdb -webfrom annovar clinvar_20240917 "${DB_DIR}"

echo ">>> 下载GWAS Catalog数据库..."
perl "${ANNOVAR_DIR}/annotate_variation.pl" -buildver "${BUILD}" -downdb -webfrom annovar gwascatalog "${DB_DIR}"

echo ">>> 下载1000 Genomes常见变异数据库..."
perl "${ANNOVAR_DIR}/annotate_variation.pl" -buildver "${BUILD}" -downdb -webfrom annovar 1000g2015aug "${DB_DIR}"

echo ">>> 下载dbSNP数据库..."
perl "${ANNOVAR_DIR}/annotate_variation.pl" -buildver "${BUILD}" -downdb -webfrom annovar avsnp150 "${DB_DIR}"

