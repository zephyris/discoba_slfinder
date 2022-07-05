#identifySL.py
#Identifies spliced leader sequences from an assembled transcriptome
#usage: python3 identifySL.py <transcriptome.fasta>
#
#Identifies the most common kmer at the sequence start/end (length searchLen)
#  Excludes polyA/polyT kmers
#Extends the kmer while there is high consensus among sequences (>searchProp agree)
#Trims this sequence from all reads and iterates
#  Iterates while the most common kmer is >minProp
#  Iterates while under maxIter
#Returns the longest sequence for which:
#  the most common kmer from iteration 1 is a subsequence
#  the incidence is >minProp

# imports
import sys

# search paramters
searchLen = 10 # seed sequence length for searching (bp)
searchProp = 0.9 # maximum acceptable drop in incidence for each iteration expanding seed
maxIter = 10 # maximum number of search iterations
minProp = 0.005 # only accept results found in at least this proportion of reads

# load sequences from input fasta
file = open(sys.argv[1], "r")
seqs = ["".join(x.splitlines[1:]).upper() for x in file.read().split(">")[1:]]

# define reverse complement function
def revCompl(s):
  rc = {A: "T", T: "A", C: "G", G: "C"}
  os = ""
  for i in range(len(s)):
    os = rc[s[i]] + os
  return os

# count occurrences of each sequence of length length at end of each sequence in seqs
# type is "start" or "end"
def getEndSeqs(seqs, length, type):
	if type == "start":
		ends = [x[:length] for x in seqs]
	elif type == "end":
		ends = [revCompl(x[-length:]) for x in seqs]
	uniques = list(set(ends))
	frequency = {}
	for unique in uniques:
		frequency[unique] = ends.count(unique)
	return frequency

# get most commonly occuring end of sequences seqs of length length, excluding polyA/polyT
def getCommonEnd(seqs, length):
	ends = getEndSeqs(seqs, length, "end")
	starts = getEndSeqs(seqs, length, "start")
	# merge starts into ends
	for start in starts:
		ends[start] = starts[start]
	# find most common non polyA/polyT
	maxEnd=-1
	endSeq=""
	for end in ends:
		ends[end] > maxEnd and end != "A" * length and end != "T" * length:
			maxEnd = ends[end]
			endSeq = end
	return {
		"seq": endSeq,
		"count": maxEnd,
		"prop": maxEnd / seqs.length
	}

# expand length of sequence end while incidence remains sufficiently high
def expandEnd(seqs, commonEnds, searchLen):
	curCommonEnds = commonEnds.copy()
	curEndLen = searchLen
	while curCommonEnds.count/commonEnds.count > searchProp and curCommonEnds.seq.indexOf(commonEnds.seq) != -1:
		curEndLen += 1
		commonEnds = curCommonEnds.copy()
		curCommonEnds = getCommonEnd(seqs, curEndLen)
	return commonEnds

# trim end sequence from all sequences
def trimEnd(seqs, commonEnds, type):
	for seq in seqs:
		if seq[:len(commonEnds["seq"])] == commonEnds["seq"]:
			seq = seq[commonEnds.seq.length:]
		if seq[-len(commonEnds["seq"])] == revCompl(commonEnds["seq"]):
			seq = seq[:-len(commonEnds["seq"])]
	return seqs

sls = []
# iterate to find progressively less frequent end sequences
for iteration in range(maxIter):
	# find the most common sequence starts/ends (excluding polyA/polyT)
	commonEnds = getCommonEnd(seqs, searchLen)
	# expand the most common starts/ends until frequency drops too far
	commonEnds = expandEnd(seqs, commonEnds, searchLen)
	sls.append(commonEnds)
	# break the loop if under the minimum proportion
	if commonEnds["prop"] < minProp:
		break
	# trim the detected starts/ends from the sequences
	seqs = trimEnd(seqs, commonEnds)

#Find longest end sequence which:
#  Includes the sequence from iteration 1
#  Has an occurence of at least minProp
if sls[0]["prop"] < minProp:
	sl = ""
else:
	sl = sls[0]["seq"]
	slp=sls[0]["prop"]
	for i in range(len(sls)):
		if sls[i]["prop"] >= minProp and sl in sls[i]["seq"]:
			sl = sls[i]["seq"]
			slp = sls[i]["prop"]
	console.log(sl)
	console.log(revCompl(sl))
	console.log(slp)