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

TRIMLINES=$(( READNUMBER*4 ))

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

# assembly
if [ -z "$2" ]; then
  # single end mode
  echo "Running in single read mode: Unpaired: $1"
  # extract reads, trim to requested number if trimlines is non-zero
  if [ "${1##*.}" == "gz" ]; then
    if [ $TRIMLINES == 0 ]; then
      gzip -cd ../$1 > U.fq
    else
      gzip -cd ../$1 | head -n $TRIMLINES > U.fq
    fi
  else
    if [ $TRIMLINES == 0 ]; then
      cat ../$1 > U.fq
    else
      cat ../$1 | head -n $TRIMLINES > U.fq
    fi
  fi
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
  Trinity --seqType fq --max_memory $MEM --single U.cor.fq --CPU $CPUS --output trinityTrial
else
  # paired end mode
  echo "Running in paired read mode: Left/Forward/1: $1, Right/Reverse/2: $2"
  # extract reads, trim to requested number if trimlines is non-zero
  if [ "${1##*.}" == "gz" ]; then
    if [ $TRIMLINES == 0 ]; then
      gzip -cd ../$1 > F.fq
    else
      gzip -cd ../$1 | head -n $TRIMLINES > F.fq
    fi
  else
    if [ $TRIMLINES == 0 ]; then
      cat ../$1 > F.fq
    else
      cat ../$1 | head -n $TRIMLINES > F.fq
    fi
  fi
  if [ "${2##*.}" == "gz" ]; then
    if [ $TRIMLINES == 0 ]; then
      gzip -cd ../$2 > R.fq
    else
      gzip -cd ../$2 | head -n $TRIMLINES > R.fq
    fi
  else
    if [ $TRIMLINES == 0 ]; then
      cat ../$2 > R.fq
    else
      cat ../$2 | head -n $TRIMLINES > R.fq
    fi
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
  #Remove spliced leader
  if [ $SPLICEDLEADERF != "N" ]; then
    cutadapt -m 25 -j $CPUS -A $SPLICEDLEADERR -G $SPLICEDLEADERF -o F.cor.sl.fq -p R.cor.sl.fq F.cor.fq R.cor.fq
    mv F.cor.sl.fq F.cor.fq
    mv R.cor.sl.fq R.cor.fq
  else
    echo "Not performing spliced leader trimming"
  fi
  #Trinity transcriptome assembly
  Trinity --seqType fq --max_memory $MEM --left F.cor.fq --right R.cor.fq --CPU $CPUS --output trinity
fi

# move files to output location and tidy
cp trinity/Trinity.fasta Trinity.fasta
cd ..
