#!/bin/bash
#SBATCH --job-name=annovar_job          # 任务名称
#SBATCH --ntasks=1                      # 总任务数
#SBATCH --cpus-per-task=4               # 每个任务分配的 CPU 核心数
#SBATCH --mem=16G                       # 分配的内存大小（如16G）
#SBATCH --output=annovar_%j.out         # 标准输出文件（%j为任务ID）
#SBATCH --error=annovar_%j.err          # 错误输出文件

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="

# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."

# 设置变量
ANNOVAR_DIR="/mnt/raid6/bacphagenetwork/tools/annovar"          # ANNVAR 的安装路径
DB_DIR="${ANNOVAR_DIR}/humandb"         # 数据库存放路径
FILE_DIR_PREFIX="/mnt/raid6/bacphagenetwork/data/13_ANNOVAR/Guangzhou"   # 文件存放路径前缀
INPUT_FILE="${FILE_DIR_PREFIX}/gz_for_annovar.csv"          # 输入文件
OUTPUT_PREFIX="${FILE_DIR_PREFIX}/results/gz_annotated_results"       # 输出文件前缀
BUILD="hg38"                            # 基因组版本

# 运行注释命令
echo ">>> Start annotating mutations..."
perl "${ANNOVAR_DIR}/table_annovar.pl" "${INPUT_FILE}" "${DB_DIR}" -buildver "${BUILD}" \
    -out "${OUTPUT_PREFIX}" -remove \
    -protocol refGene,clinvar_20240917,avsnp150,ALL.sites.2015_08 \
    -operation g,f,f,f -nastring . -csvout

echo ">>> Annotation completes!The result is saved at ${OUTPUT_PREFIX}.${BUILD}_multianno.csv"
