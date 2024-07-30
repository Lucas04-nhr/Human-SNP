#! /bin/bash

# Check if the input folder exists

echo "Checking the input folder..."

if [ ! -d "/mnt/raid6/bacphagenetwork/data/03_samtools_sorted" ]
then
    echo "The input folder does not exist, please check manually ..."
    return 1
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

mkdir log
mkdir log/Beijing
mkdir log/Guangzhou
echo "The log folder has been created."

# Create the output folder

echo "Creating the output folder..."

if [ -d "/mnt/raid6/bacphagenetwork/data/04_samtools_tag_added" ]
then
    echo "The output folder already exists, do you want to back up the output folder? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Backing up the output folder ..."
        mv /mnt/raid6/bacphagenetwork/data/03_samtools_collate /mnt/raid6/bacphagenetwork/data/03_samtools_collate_bak_$(date +%Y%m%d%H%M%S)
    else
        echo "The original output folder will be removed."
        rm -rf /mnt/raid6/bacphagenetwork/data/04_samtools_tag_added
    fi
fi

mkdir /mnt/raid6/bacphagenetwork/data/04_samtools_tag_added
mkdir /mnt/raid6/bacphagenetwork/data/04_samtools_tag_added/Beijing
mkdir /mnt/raid6/bacphagenetwork/data/04_samtools_tag_added/Guangzhou

echo "The output folder has been created."

# Check if the requirements satisfied

echo "Checking the requirements ..."

if ! command -v samtools &> /dev/null
then
    echo "samtools is not installed. Please install samtools and try again."
    return 1
else
    echo "samtools is installed."
fi

if ! command -v bedtools &> /dev/null
then
    echo "bedtools is not installed. Please install bedtools and try again."
    return 1
else
    echo "bedtools is installed."
fi

echo "All requirements are satisfied."

# Creating sbatch files

echo "Creating sbatch files ..."

bash ./create_sbatch_BJ.sh
bash ./create_sbatch_GZ.sh

echo "The sbatch files have been created."