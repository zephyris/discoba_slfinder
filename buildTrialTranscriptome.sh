#!/bin/bash

# builds a trial transcriptome from (optionally gz compressed) fastqs of illumina reads
# accepts a single file for either unpaired or two files, one for each paired direction
# for paired reads: buildTrialTranscriptome.sh forward.fastq[.gz] reverse.fastq[.gz]
# for single reads: buildTrialTranscriptome.sh reads.fastq[.gz]

# global variable defining CPU and memory availability
# automatic all CPU threads
CPUS=$(nproc --all)
# automatic all free RAM
MEM=$(free --giga|awk '/^Mem:/{print $4}')G
echo "$CPUS CPUs, $MEM RAM"

# get script directory
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# paths for assembly packages
RCORRECTOR=$SCRIPTDIR/rcorrector/run_rcorrector.pl
FILTERUNCORR=$SCRIPTDIR/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py
TRIMGALORE=$SCRIPTDIR/trimgalore-v0.6.0/trim_galore
export PATH=$PATH:$SCRIPTDIR/salmon-v1.5.2/bin
export PATH=$PATH:$SCRIPTDIR/trinityrnaseq-v2.12.0

# check if assembly packages seem to be installed
# nb. this doesn't check if each can actually run1
if [ ! -f $SCRIPTDIR/versions_assembly.txt ]; then
  echo "Is the transcriptome build software correctly installed?"
  echo "Please run installTranscriptomeBuilding.sh"
  exit
fi

echo "Building trial transcriptome using program versions:"
cat $SCRIPTDIR/versions_assembly.txt

if [ -z "$1" ] && [ -z "$2" ]; then
  echo "No FASTQ files given"
  exit
fi

# assembly of a subsample of 1 million reads
# sl can be identified from relatively few mRNAs
if [ -z "$2" ]; then
  # single end mode
  echo "Running in single read mode: Unpaired: $1"
  if [ "${1##*.}" == "gz" ]; then
    gzip -cd $1 | head -n 40000000 > U.fq
  else
    cat $1 | head -n 40000000 > U.fq
  fi
  #Use standard correction, filtering, trimming
  perl $RCORRECTOR -t $CPUS -s U.fq
  $TRIMGALORE --output_dir tgl --length 36 -q 5 --stringency 1 -e 0.1 U.cor.fq
  mv tgl/U.cor_trimmed.fq U.cor.fq
  rm -r tgl
  #Trial Trinity transcriptome assembly
  Trinity --seqType fq --max_memory $MEM --single U.cor.fq --CPU $CPUS --output trinityTrial
else
  # paired end mode
  echo "Running in paired read mode: Forward/1: $1, Reverse/2: $2"
  if [ "${1##*.}" == "gz" ]; then
    gzip -cd $1 | head -n 40000000 > F.fq
  else
    cat $1 | head -n 40000000 > F.fq
  fi
  if [ "${2##*.}" == "gz" ]; then
    gzip -cd $2 | head -n 40000000 > R.fq
  else
    cat $2 | head -n 40000000 > R.fq
  fi
  #Use standard correction, filtering, trimming
  perl $RCORRECTOR -t $CPUS -1 F.fq -2 R.fq
  python $FILTERUNCORR -1 F.cor.fq -2 R.cor.fq -s cor
  $TRIMGALORE --paired --retain_unpaired --phred33 --output_dir tgl --length 36 -q 5 --stringency 1 -e 0.1 unfixrm_F.cor.fq unfixrm_R.cor.fq
  mv tgl/unfixrm_F.cor_val_1.fq F.cor.fq
  mv tgl/unfixrm_R.cor_val_2.fq R.cor.fq
  mv tgl/unfixrm_F.cor_unpaired_1.fq U.cor.fq
  cat tgl/unfixrm_R.cor_unpaired_2.fq >> U.cor.fq
  rm tgl/unfixrm_R.cor_unpaired_2.fq
  rm -r tgl
  #Trial Trinity transcriptome assembly
  Trinity --seqType fq --max_memory $MEM --left F.cor.fq --right R.cor.fq --CPU $CPUS --output trinityTrial
fi
