#! /bin/bash

ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis'

cd $ANALYSIS_PATH

echo "Current directory: $(pwd)"

echo "Moving files..."

if (! [ -d ./Beijing ])
then
    mkdir ./Beijing
fi

if (! [ -d ./Guangzhou ])
then
    mkdir ./Guangzhou
fi

mv ./BJ*.sam ./Beijing
mv ./GZ*.sam ./Guangzhou

echo "*.sam files have been moved."

cd log

echo "Current directory: $(pwd)"

echo "Moving files..."

if (! [ -d ./Beijing ])
then
    mkdir ./Beijing
fi

if (! [ -d ./Guangzhou ])
then
    mkdir ./Guangzhou
fi

mv ./BJ*.log ./Beijing
mv ./GZ*.log ./Guangzhou

echo "*.log files have been moved."

echo "Great! All files have been moved."
