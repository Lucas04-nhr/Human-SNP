#! /bin/bash

# Check if the input folder exists

echo "Checking the input folder..."

if [ ! -d "/mnt/raid6/bacphagenetwork/data/04_dulplicate_marked" ]
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

if [ -d "/mnt/raid6/bacphagenetwork/data/05_breakpoint_marked" ]
then
    echo "The output folder already exists, do you want to back up the output folder? (y/n)"
    read answer
    if [ $answer == "y" ] || [ $answer == "Y" ]
    then
        echo "Backing up the output folder ..."
        mv /mnt/raid6/bacphagenetwork/data/05_breakpoint_marked /mnt/raid6/bacphagenetwork/data/05_breakpoint_marked_bak_$(date +%Y%m%d%H%M%S)
    else
        echo "The original output folder will be removed."
        rm -rf /mnt/raid6/bacphagenetwork/data/05_breakpoint_marked
    fi
fi

mkdir /mnt/raid6/bacphagenetwork/data/05_breakpoint_marked
mkdir /mnt/raid6/bacphagenetwork/data/05_breakpoint_marked/Beijing
mkdir /mnt/raid6/bacphagenetwork/data/05_breakpoint_marked/Guangzhou

echo "The output folder has been created."

# Check if the requirements satisfied

echo "Checking the requirements ..."

# System-wide java
if [ ! -f "/bin/java" ]
then
    echo "The system-wide java is not found, please check manually ..."
    exit 1
fi

# JDK22.0.1
if [ ! -f "/mnt/raid6/bacphagenetwork/tools/jdk-22.0.1/bin/java" ]
then
    echo "The external JDK is not found, please check manually ..."
    exit 1
fi

# GATK 4.3
if [ ! -f "/mnt/raid6/bacphagenetwork/tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar" ]
then
    echo "The GATK v4.3 is not found, please check manually ..."
    exit 1
fi

# GATK4.5
if [ ! -f "/mnt/raid6/bacphagenetwork/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" ]
then
    echo "The GATK v4.5 is not found, please check manually ..."
    exit 1
fi

# Picard 2.26
if [ ! -f "/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v2.26.0/picard.jar" ]
then
    echo "The Picard v2.26 is not found, please check manually ..."
    exit 1
fi

# Picard 3.0
if [ ! -f "/mnt/raid6/bacphagenetwork/tools/picard_pre-built/v3.0/picard.jar" ]
then
    echo "The Picard v3.0 is not found, please check manually ..."
    exit 1
fi

echo "All requirements are satisfied."

# Ask for creating sbatch files

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