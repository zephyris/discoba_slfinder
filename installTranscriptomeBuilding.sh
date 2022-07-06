#!/bin/bash

if [ "$EUID" -ne 0 ]; then
  echo "Please run as root"
  exit
fi

# apt-get installations
sudo apt-get update
# general packages
sudo apt-get --yes install curl python3 python3-pip git
# assembly packages
sudo apt-get --yes install cutadapt jellyfish bowtie2 samtools
# dependencies for building trinity
sudo apt-get --yes install autoconf libbz2-dev liblzma-dev

# git Trinity and compile
if [ ! -d trinityrnaseq-v2.12.0 ]; then
  curl https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.12.0/trinityrnaseq-v2.12.0.FULL.tar.gz -L -o trinityrnaseq.tar.gz
  tar -xvf trinityrnaseq.tar.gz
  cd trinityrnaseq-v2.12.0
  make
  cd ..
  rm trinityrnaseq.tar.gz
fi

# download Salmon and extract
if [ ! -d salmon-v1.5.2 ]; then
  curl https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz -L -o salmon.tar.gz
  tar -xvf salmon.tar.gz
  mv salmon-1.5.2_linux_x86_64 salmon-v1.5.2
  rm salmon.tar.gz
fi

# git RCorrector and compile
if [ ! -d rcorrector ]; then
  git clone https://github.com/mourisl/rcorrector.git
  cd rcorrector
  make
  cd ..
fi

# git TranscriptomeAssemblyTools
if [ ! -d TranscriptomeAssemblyTools ]; then
  git clone https://github.com/harvardinformatics/TranscriptomeAssemblyTools
fi

# download TrimGalore! and extract
if [ ! -d trimgalore-v0.6.0 ]; then
  curl https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.zip -L -o trimgalore.zip
  unzip trimgalore.zip
  mv TrimGalore-0.6.0 trimgalore-v0.6.0
  rm trimgalore.zip
fi

# Paths
SCRIPTDIR=$(pwd)
RCORRECTOR=$SCRIPTDIR/rcorrector/run_rcorrector.pl
FILTERUNCORR=$SCRIPTDIR/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py
TRIMGALORE=$SCRIPTDIR/trimgalore-v0.6.0/trim_galore
export PATH=$PATH:$SCRIPTDIR/salmon-v1.5.2/bin
export PATH=$PATH:$SCRIPTDIR/trinityrnaseq-v2.12.0

# Check everything runs and record versions
# If versions_assembly.txt is lacking a version for anything then it hasn't installed correctly
echo ""
echo "Installed versions:"
echo "rcorrector" $(perl $RCORRECTOR -version | sed "s/v//" | awk '{print $2}') > versions_assembly.txt
echo "trimgalore" $($TRIMGALORE -v | grep "version" | sed "s/version//g") >> versions_assembly.txt
echo "trinity" $(Trinity --version | grep "Trinity version:" | sed "s/Trinity version: Trinity-v//g") >> versions_assembly.txt
salmon --version >> versions_assembly.txt
jellyfish --version >> versions_assembly.txt
echo "cutadapt" $(cutadapt --version) >> versions_assembly.txt
echo "bowtie2" $(bowtie2 --version | head -n 1 | awk '{print $3}') >> versions_assembly.txt
cat versions_assembly.txt

