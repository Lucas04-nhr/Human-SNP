#! /bin/bash

# Check if the source data folder exists

export DATA_PATH=/mnt/raid6/bacphagenetwork/data/
echo "Checking the source data folder..."

if [ ! -d "${DATA_PATH}/skin_metagenome/" ]
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

mkdir -p log/01_align/Beijing
mkdir -p log/01_align/Guangzhou
mkdir -p log/02_filter/Beijing
mkdir -p log/02_filter/Guangzhou
mkdir -p log/03_sort/Beijing
mkdir -p log/03_sort/Guangzhou
mkdir -p log/04_convert/Beijing
mkdir -p log/04_convert/Guangzhou
echo "The log folder has been created."

# Create the output folder

echo "Creating the output folder..."

mkdir -p ${DATA_PATH}/01_align/Beijing
mkdir -p ${DATA_PATH}/01_align/Guangzhou
mkdir -p ${DATA_PATH}/02_filter/Beijing
mkdir -p ${DATA_PATH}/02_filter/Guangzhou
mkdir -p ${DATA_PATH}/03_sort/Beijing
mkdir -p ${DATA_PATH}/03_sort/Guangzhou
mkdir -p ${DATA_PATH}/04_convert/Beijing
mkdir -p ${DATA_PATH}/04_convert/Guangzhou

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
