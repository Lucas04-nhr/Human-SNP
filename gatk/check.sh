#! /bin/bash

# Check if the source data folder exists

export DATA_PATH=/mnt/raid6/bacphagenetwork/data/
echo "Checking the source data folder..."

if [ ! -d "${DATA_PATH}/04_convert/" ]
then
    echo "The input folder does not exist, please check manually ..."
    exit 1
else
    echo "The input folder already exists."
fi

# Create the log folder

echo "Creating the log folder..."

if [ -d "./log" ]
then
    echo "The log folder already exists, do you want to back up the log folder? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Backing up the log folder ..."
        mv log log_bak_$(date +%Y%m%d%H%M%S)
    else
        echo "The original log folder will be removed."
        rm -rf log
    fi
fi

mkdir -p log/Beijing
mkdir -p log/Guangzhou
echo "The log folder has been created."

# Create the output folder

echo "Creating the output folder..."

mkdir -p ${DATA_PATH}/05_BaseRecalibrator/Beijing
mkdir -p ${DATA_PATH}/05_BaseRecalibrator/Guangzhou
mkdir -p ${DATA_PATH}/06_ApplyBQSR/Beijing
mkdir -p ${DATA_PATH}/06_ApplyBQSR/Guangzhou
mkdir -p ${DATA_PATH}/07_HaplotypeCaller/Beijing
mkdir -p ${DATA_PATH}/07_HaplotypeCaller/Guangzhou


echo "Do you want to create sbatch files? (y/n)"
read answer
if [ $answer != "y" ] && [ $answer != "Y" ]
then
    echo "The process has been terminated."
    exit 0
fi

# Creating sbatch files

echo "Creating sbatch files ..."

bash ./create_sbatch_BJ.sh
bash ./create_sbatch_GZ.sh

echo "The sbatch files have been created."
