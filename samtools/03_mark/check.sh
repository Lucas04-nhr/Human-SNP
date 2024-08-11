#! /bin/bash

# Check if the input folder exists

echo "Checking the input folder..."

if [ ! -d "/mnt/raid6/bacphagenetwork/data/03_samtools_sorted" ]
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

mkdir log
mkdir log/Beijing
mkdir log/Guangzhou
echo "The log folder has been created."

# Create the output folder

echo "Creating the output folder..."

if [ -d "/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked" ]
then
    echo "The output folder already exists, do you want to back up the output folder? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Backing up the output folder ..."
        mv /mnt/raid6/bacphagenetwork/data/04_dulplicate_marked /mnt/raid6/bacphagenetwork/data/04_dulplicate_marked_bak_$(date +%Y%m%d%H%M%S)
    else
        echo "The original output folder will be removed."
        rm -rf /mnt/raid6/bacphagenetwork/data/04_dulplicate_marked
    fi
fi

mkdir /mnt/raid6/bacphagenetwork/data/04_dulplicate_marked
mkdir /mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Beijing
mkdir /mnt/raid6/bacphagenetwork/data/04_dulplicate_marked/Guangzhou

echo "The output folder has been created."

# Check if the requirements satisfied

echo "Checking the requirements ..."

export JAVA_BIN='/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java'
export PICARD_BIN='/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar'

if [ ! -f "$JAVA_BIN" ]
then
    echo "Error: $JAVA_BIN does not exist."
    exit 1
else
    echo "The file $JAVA_BIN exists."
fi

if [ ! -f "$PICARD_BIN" ]
then
    echo "Error: $PICARD_BIN does not exist."
    exit 1
else
    echo "The file $PICARD_BIN exists."
fi

echo "All requirements are satisfied."

# Creating sbatch files

echo "Creating sbatch files ..."

bash ./create_sbatch_BJ.sh
bash ./create_sbatch_GZ.sh

echo "The sbatch files have been created."