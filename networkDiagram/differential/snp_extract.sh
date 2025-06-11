#! /bin/bash
#SBATCH --job-name=snp_extract
#SBATCH --output=./snp_extract.%j.out
#SBATCH --error=./snp_extract.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="

# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="

#plink可执行路径
export PATH="/mnt/raid6/bacphagenetwork/tools/plink_1.9_linux_x86_64:$PATH"

# 设置变量
PLINK_PREFIX="/mnt/raid6/bacphagenetwork/data/12_plink/Full/output/bac_age/converted_genotyped"
SNP_LIST="/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/networkDiagram/differential/snp_list.txt"
OUT_PREFIX="/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/networkDiagram/differential/output/snp_extract"
INPUT_DIR="/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/networkDiagram/differential/output"

mkdir -p "$INPUT_DIR"

# 写入目标 SNP 到 snplist
echo "Extracting SNPs listed in $SNP_LIST..."
plink \
  --bfile "$PLINK_PREFIX" \
  --extract "$SNP_LIST" \
  --recodeA \
  --allow-extra-chr \
  --out "$OUT_PREFIX"

echo "目标SNP已提取，结果保存在 $OUT_PREFIX.raw"
