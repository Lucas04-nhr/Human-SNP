#!/bin/bash

# 读取 BJ_sbatch.list 文件中的每一行
while IFS= read -r infile; do
    # 提取 sample_name
    sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

    # 设置 MANTA_ANALYSIS_PATH 环境变量
    export PARENT_PATH='/mnt/raid6/bacphagenetwork/data/07_manta/01_exec/Beijing'
    export MANTA_ANALYSIS_PATH="$PARENT_PATH/$sample_name"
    echo "The path to the execute folder of $sample_name has been set to $MANTA_ANALYSIS_PATH."

    echo "Initializing complete."
    echo "========================================"

    # 复制 Manta 分析结果
    echo "Copying the Manta analysis result of ${sample_name}..."
    export STATS_PATH="${MANTA_ANALYSIS_PATH}/results"
    export RESULT_PATH="/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/manta/02_result/Beijing"

    mkdir -p "${RESULT_PATH}"

    if [ -d "${RESULT_PATH}/${sample_name}" ]; then
        echo "Error: The result folder already exists, overwriting it."
        rm -rf "${RESULT_PATH}/${sample_name}" 
    fi
    
    cp -r "${STATS_PATH}" "${RESULT_PATH}" \
    || { echo "Error: Copying the Manta analysis result failed."; exit 1; }
    echo "Copying complete."

    echo "Renaming the result folder..."
    mv "${RESULT_PATH}/results" "${RESULT_PATH}/${sample_name}" \
    || { echo "Error: Renaming the result folder failed."; exit 1; }
    echo "Renaming complete."

    echo "========================================"
    echo "The Manta analysis result of ${sample_name} has been successfully copied to ${RESULT_PATH}/${sample_name}."

done < BJ_sbatch.list

echo "All Manta analysis results have been successfully copied."