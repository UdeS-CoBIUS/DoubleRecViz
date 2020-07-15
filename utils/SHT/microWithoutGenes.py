#!/usr/bin/python

'''
Author : Marie Degen

Labeled all the exons with his corresponding alternative splicing event from a gene family. 
Example : python microWithoutGenes.py FAM86_microalignment.fasta

'''

import sys
import argparse
import os
from Bio import SeqIO

#Save Path 

def microWithoutGenes(microalignmentFile):
	path =  sys.path[0]	
	save_path = path + "/microalignmentWithoutGene/"	

	tmp = microalignmentFile.split("/")
	tmp = tmp[-1]
	tmp = tmp.split("_")[0]
	source2target = open("initialSource/" +  tmp + "_initialsource2target.txt", "r")
	lines = source2target.readlines()
	transcripts = []
	genes =[]
	for line in lines:
		line = line.replace("\n", "")
		parts = line.split(" ")
		
		transcripts.append(parts[0])
		genes.append(parts[1])


	dictSequence = {}

	for record in SeqIO.parse(microalignmentFile, "fasta"):
		dictSequence[record.id] = str(record.seq)

	namefile = microalignmentFile.split("/")[-1]

	completeName = os.path.join(save_path, namefile)         
	outputfile = open(completeName, "w")

	for k, v in dictSequence.items():
		
		if k in transcripts:
			outputfile.write(">" + k + "\n" + v + "\n")
	outputfile.close()

	return completeName


