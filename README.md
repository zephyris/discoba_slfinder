# discoba_slfinder
Simple tools for automated identification of spliced leader sequences from transcriptome sequencing.

slfinder is a simple, and quite naive, tool for finding the most common sequence found in whole or in part at the start or end (reverse complemented) of sequences in a fasta file, excluding poly A. When run on a transcriptome in which an extremely large proportion of transcripts have a spliced leader from a common mini-exon this should identify the spliced leader sequence. This simple implementation is particularly suited to the many Discoba species which use a single spliced leader for all mRNAs.
The reported sequence does not remove any nucleotides which always follow the splice site and do not come from the mini-exon, expected to be `TG` in kinetoplastids.

## Usage
Two implementations of slfinder are included, one written in python, one in nodejs[^1].
[^1]: The nodejs implementation is easily an order of magnitude faster than python3 and gives an identical output.


Run using either:

`python3 ../path/to/discoba_slfinder/slfinder.py transcriptome.fasta` or `nodejs ../path/to/discoba_slfinder/slfinder.js transcriptome.fasta`

A verbose output is displayed in the console, for example:
```
Seed kmer length: 10bp
Search consensus proportion: 0.9
Maximum iterations: 10
Minimum proportion: 0.005
Number of sequences: 14578
iter.   type            length  count   perc.   sequence
0       seed seq.       10bp    2392    16.41%  GAACAGTTTC
0       expanded        22bp    2392    16.41%  GAACAGTTTCTGTACTATATTG
1       seed seq.       10bp    1295    8.88%   AGAACAGTTT
1       expanded        23bp    1295    8.88%   AGAACAGTTTCTGTACTATATTG
2       seed seq.       10bp    1236    8.48%   AACAGTTTCT
2       expanded        21bp    1236    8.48%   AACAGTTTCTGTACTATATTG
3       seed seq.       10bp    772     5.3%    TAGAACAGTT
3       expanded        24bp    770     5.28%   TAGAACAGTTTCTGTACTATATTG
4       seed seq.       10bp    679     4.66%   AGTTTCTGTA
4       expanded        18bp    678     4.65%   AGTTTCTGTACTATATTG
5       seed seq.       10bp    641     4.4%    AATAGAACAG
5       expanded        26bp    641     4.4%    AATAGAACAGTTTCTGTACTATATTG
6       seed seq.       10bp    549     3.77%   CAGTTTCTGT
6       expanded        19bp    548     3.76%   CAGTTTCTGTACTATATTG
7       seed seq.       10bp    481     3.3%    GTTTCTGTAC
7       expanded        17bp    481     3.3%    GTTTCTGTACTATATTG
8       seed seq.       10bp    230     1.58%   ACAGTTTCTG
8       expanded        20bp    230     1.58%   ACAGTTTCTGTACTATATTG
9       seed seq.       10bp    217     1.49%   ATAGAACAGT
9       expanded        25bp    217     1.49%   ATAGAACAGTTTCTGTACTATATTG
Updated from iteration 2 AGAACAGTTTCTGTACTATATTG
Updated from iteration 4 TAGAACAGTTTCTGTACTATATTG
Updated from iteration 6 AATAGAACAGTTTCTGTACTATATTG
Spliced leader sequence identified:
AATAGAACAGTTTCTGTACTATATTG
Reverse complement:
CAATATAGTACAGAAACTGTTCTATT
Frequency of full length sequence:
0.04397
Frequency of whole or partial sequence:
0.58225
```

To record the output, pipe the output to a file. For example:
`python3 slfinder.py transcriptome.fasta > sl.txt`

Only a minimal output of the identified spliced leader sequence, its reverse complement and its whole and whole or in part frequencies will be written to the file. In this example, `sl.txt` will contain:
```
AATAGAACAGTTTCTGTACTATATTG
CAATATAGTACAGAAACTGTTCTATT
0.04397
0.58225
```

## Transcriptome
A fasta file of a de-novo assembled transcriptome is required as an input. A relatively low-quality transcriptome, ie. only the more abundant mRNAs, is OK and can be generated relatively quickly from a small number (eg. 1 million) of Illumina RNAseq reads. In general, transcriptomes built from longer reads and paired ends are better for spliced leader detection.

`buildTrialTranscriptome.sh` is provided to perform a simple de-novo build from Illumina `.fastq` or `.fastq.gz` files.
To use, first install the necessary dependencies:

`sudo bash installTranscriptomeBuilding.sh`

This will install necessary packages from apt-get (cutadapt, jellyfish, bowtie2, samtools) and download and build or extract transcript and transcriptome tools (trinity, salmon, rcorrector, transcriptomeassemblytools and trimgalore). This places the built/extracted tools in the same directory as `installTranscriptomeBuilding.sh` and `buildTrialTranscriptome.sh`.

To build a transcriptome run:

`bash ../path/to/discoba_slfinder/buildTrialTranscriptome.sh pairedreads_1.fastq.gz pairedreads_2.fastq.gz` or `bash ../path/to/discoba_slfinder/buildTrialTranscriptome.sh unpairedreads.fastq`

The trial transcriptome is output to `trinityTrial/Trinity.fasta`.
