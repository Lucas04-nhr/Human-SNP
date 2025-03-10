#!/bin/bash
#SBATCH --job-name=plink_filter_job  # 作业名称
#SBATCH --output=plink_filter_output.log  # 标准输出和错误日志文件
#SBATCH --mem=4G  # 内存需求

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="
echo ""

# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="

# 运行Python脚本
python plinkFilter.py data/12_plink/Full/modified/bac_age Full/merged/filtered_bac_age.csv
echo "The plinkFilter.py has been executed."