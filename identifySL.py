#!/usr/bin/python3

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

show_progress = True

def print_status(s):
  if show_progress == True:
    sys.stderr.write(str(s)+"\n")

# search paramters
searchLen = 10 # seed sequence length for searching (bp)
maxSearchLen = 40 # maximum seed sequence length to search for (bp)
searchProp = 0.9 # maximum acceptable drop in incidence for each iteration expanding seed
maxIter = 10 # maximum number of search iterations
minProp = 0.005 # only accept results found in at least this proportion of reads

print_status("Seed kmer length: "+str(searchLen)+"bp");
print_status("Max kmer length: "+str(maxSearchLen)+"bp");
print_status("Search consensus proportion: "+str(searchProp));
print_status("Maximum iterations: "+str(maxIter));
print_status("Minimum proportion: "+str(minProp));

# load sequences from input fasta
file = open(sys.argv[1], "r")
seqs = ["".join(x.splitlines()[1:]).upper() for x in file.read().split(">")[1:]]

# define reverse complement function
def revCompl(s):
  rc = {"A": "T", "T": "A", "C": "G", "G": "C"}
  os = ""
  for i in range(len(s)):
    os = rc[s[i]] + os
  return os

# count occurrences of each sequence of length length at end of each sequence in seqs
# type is "start" or "end"
def getEndSeqs(seqs, length, type):
	# only analyse sequences longer than length
	filteredSeqs = filter(lambda x: len(x) >= length, seqs)
	if type == "start":
		ends = [x[:length] for x in filteredSeqs]
	elif type == "end":
		ends = [revCompl(x[-length:]) for x in filteredSeqs]
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
		if start in ends:
			ends[start] += starts[start]
		else:
			ends[start] = starts[start]
	# find most common non polyA/polyT
	maxEnd=-1
	endSeq=""
	polyA = "A" * length
	polyT = "T" * length
	for end in ends:
		if ends[end] > maxEnd and end != polyA and end != polyT:
			maxEnd = ends[end]
			endSeq = end
	return {
		"seq": endSeq,
		"count": maxEnd,
		"prop": maxEnd / float(len(seqs))
	}

# expand length of sequence end while incidence remains sufficiently high
def expandEnd(seqs, commonEnds, searchLen):
	curCommonEnds = commonEnds.copy()
	curEndLen = searchLen
	while curCommonEnds["count"] / commonEnds["count"] > searchProp and commonEnds["seq"] in curCommonEnds["seq"] and curEndLen <= maxSearchLen:
		curEndLen += 1
		commonEnds = curCommonEnds.copy()
		curCommonEnds = getCommonEnd(seqs, curEndLen)
	return commonEnds

# trim end sequence from all sequences
def trimEnd(seqs, commonEnds):
	for s in range(len(seqs)):
		if seqs[s][:len(commonEnds["seq"])] == commonEnds["seq"]:
			seqs[s] = seqs[s][len(commonEnds["seq"]):]
		if seqs[s][-len(commonEnds["seq"]):] == revCompl(commonEnds["seq"]):
			seqs[s] = seqs[s][:-len(commonEnds["seq"])]
	return seqs

sls = []
# iterate to find progressively less frequent end sequences
print_status("\t".join(["iter.", "type    ", "length", "count", "perc.", "sequence"]));
seqsWorking = seqs.copy()
for iteration in range(maxIter):
	# find the most common sequence starts/ends (excluding polyA/polyT)
	commonEnds = getCommonEnd(seqsWorking, searchLen)
	print_status("\t".join([str(x) for x in [iteration, "seed seq.", str(len(commonEnds["seq"]))+"bp", commonEnds["count"], str(round(100*commonEnds["prop"], 2))+"%", commonEnds["seq"]]]))
	# expand the most common starts/ends until frequency drops too far
	commonEnds = expandEnd(seqsWorking, commonEnds, searchLen)
	print_status("\t".join([str(x) for x in [iteration, "expanded", str(len(commonEnds["seq"]))+"bp", commonEnds["count"], str(round(100*commonEnds["prop"], 2))+"%", commonEnds["seq"]]]))
	sls.append(commonEnds)
	# break the loop if under the minimum proportion
	if commonEnds["prop"] < minProp:
		break
	# trim the detected starts/ends from the sequences
	seqsNew = trimEnd(seqsWorking, commonEnds)
	seqsWorking = seqsNew.copy()

#Find longest end sequence which:
#  Includes the sequence from iteration 1
#  Has an occurence of at least minProp
if sls[0]["prop"] < minProp:
	sl = ""
	print_status("No common start/end with necessary occurence found")
else:
	sl = sls[0]["seq"]
	slp=sls[0]["prop"]
	for i in range(len(sls)):
		if sls[i]["prop"] >= minProp and sl in sls[i]["seq"]:
			sl = sls[i]["seq"]
			slp = sls[i]["prop"]
	print_status("Spliced leader sequence identified:")
	print(sl)
	print_status("Reverse complement:")
	print(revCompl(sl))
	print_status("Frequency:")
	print(str(round(slp, 5)))
	cumulativeProp = 0
	for csl in sls:
		if csl["seq"] in sl:
			cumulativeProp += csl["prop"]
	print_status("Found on "+str(round(100*cumulativeProp, 2))+"% of transcripts")
