#!/bin/bash

# builds a transcriptome from (optionally gz compressed) fastqs of illumina reads
# accepts a single file for either unpaired or two files, one for each paired direction
# for paired reads: buildTranscriptome.sh forward.fastq[.gz] reverse.fastq[.gz]
# for single reads: buildTranscriptome.sh reads.fastq[.gz]
# use -n to subsample reads to this number (default all reads)
# use -s to set a spliced leader to trim from input reads (default no trimming)
# use -o to set output directory (default trinity)
# use -c to set number of CPUs (default: all CPUs)
# use -m to set RAM to use in GB (default: all free RAM)

# defaults
# automatic all CPU threads
CPUS=$(nproc --all)
# automatic all free RAM
MEM=$(free --giga|awk '/^Mem:/{print $4}')G
# settings
READNUMBER=0
SPLICEDLEADERF="N"
OUTDIR="trinity"

# variable for remaining positional arguments
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -n|--numreads)
      READNUMBER="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--splicedleader)
      SPLICEDLEADERF="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift # past argument
      shift # past value
      ;;
    -c|--cpus)
      CPUS="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--memgb)
      MEM="$2G"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

echo "$CPUS CPUs, $MEM RAM"
echo "$READNUMBER reads"

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

echo "Building transcriptome using program versions:"
cat $SCRIPTDIR/versions_assembly.txt

SPLICEDLEADERF=${SPLICEDLEADERF^^}
SPLICEDLEADERR=$(echo $SPLICEDLEADERF | rev | tr ATCG TAGC)
echo "Spliced leader:" $SPLICEDLEADERF
echo "Reverse complement:" $SPLICEDLEADERR

if [ -z "$1" ] && [ -z "$2" ]; then
  echo "No FASTQ files given"
  exit
fi

mkdir $OUTDIR
cd $OUTDIR

# samplereads input.fq(.gz) number output.fa
samplereads() {
  echo "Sampling $2 reads from $1, to $3"
  if [ $2 == 0 ]; then
    echo "No sampling, use all reads"
    if [ "${1##*.}" == "gz" ]; then
      gzip -d -c ../$1 > $3
    else
      cat ../$1 > $3
    fi
  else
    # ?just use seqtk sample two pass mode?
    # number of reads available
    echo "Counting reads available..."
    if [ "${1##*.}" == "gz" ]; then
      READCOUNT=$(gzip -d -c ../$1 | wc -l)
    else
      READCOUNT=$(cat ../$1 | wc -l)
    fi
    echo "$READCOUNT reads available"
    # if subsampling is required
    if [ $READCOUNT -gt $2 ]; then
      # calculate subsampling factor (_much_ lower memory use than set read number)
      SAMPLEFAC=$(echo "scale=8; ($2/$READCOUNT)" | bc | sed 's/^\./0./')
      echo "Subsample factor: $SAMPLEFAC"
      # do subsampling with fixed seed
      if [ "${1##*.}" == "gz" ]; then
        gzip -d -c ../$1 | seqtk sample -s1337 - $2 > $3
      else
        cat ../$1 | seqtk sample -s1337 - $2 > $3
      fi
    else
      # otherwise, just output everything
      echo "More reads requested than available, using all reads"
      if [ "${1##*.}" == "gz" ]; then
        gzip -d -c ../$1 > $3
      else
        cat ../$1 > $3
      fi
    fi
  fi
}

# samplereads input.fq(.gz) number output.fa
samplereads_simple() {
  echo "Sampling $2 reads from $1, to $3"
  if [ $2 == 0 ]; then
    echo "No sampling, use all reads"
    if [ "${1##*.}" == "gz" ]; then
      gzip -d -c ../$1 > $3
    else
      cat ../$1 > $3
    fi
  else
    echo "Using first $2 reads"
    READNUMBER=$2
    TRIMLINES=$(( READNUMBER*4 ))
    if [ "${1##*.}" == "gz" ]; then
      gzip -cd ../$1 | head -n $TRIMLINES > U.fq
    else
      cat ../$1 | head -n $TRIMLINES > U.fq
    fi
  fi
}

# assembly
if [ -z "$2" ]; then
  # single end mode
  echo "Running in single read mode: Unpaired: $1"
  # extract reads, trim to requested number if readnumber is non-zero
  samplereads $1 $READNUMBER U.fq
  #Use standard correction, filtering, trimming
  perl $RCORRECTOR -t $CPUS -s U.fq
  $TRIMGALORE --output_dir tgl --length 36 -q 5 --stringency 1 -e 0.1 U.cor.fq
  mv tgl/U.cor_trimmed.fq U.cor.fq
  rm -r tgl
  #Remove spliced leader
  if [ $SPLICEDLEADERF != "N" ]; then
    cutadapt -m 25 -j $CPUS -A $SPLICEDLEADERR -G $SPLICEDLEADERF -o U.cor.sl.fq U.cor.fq
    mv U.cor.sl.fq U.cor.fq
  else
    echo "Not performing spliced leader trimming"
  fi
  #Trinity transcriptome assembly
  Trinity --seqType fq --max_memory $MEM --single U.cor.fq --CPU $CPUS --output $OUTDIR --jaccard_clip
else
  # paired end mode
  echo "Running in paired read mode: Left/Forward/1: $1, Right/Reverse/2: $2"
  # extract reads, trim to requested number if readnumber is non-zero
  samplereads $1 $READNUMBER F.fq
  samplereads $2 $READNUMBER R.fq
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
  #Remove spliced leader
  if [ $SPLICEDLEADERF != "N" ]; then
    cutadapt -m 25 -j $CPUS -A $SPLICEDLEADERR -G $SPLICEDLEADERF -o F.cor.sl.fq -p R.cor.sl.fq F.cor.fq R.cor.fq
    mv F.cor.sl.fq F.cor.fq
    mv R.cor.sl.fq R.cor.fq
  else
    echo "Not performing spliced leader trimming"
  fi
  #Trinity transcriptome assembly
  Trinity --seqType fq --max_memory $MEM --left F.cor.fq --right R.cor.fq --CPU $CPUS --output $OUTDIR --jaccard_clip
fi

# move files to output location and tidy
cp trinity/Trinity.fasta Trinity.fasta
cd ..
