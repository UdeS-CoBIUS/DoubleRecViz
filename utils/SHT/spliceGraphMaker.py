#!/usr/bin/python

'''
Author : Marie Degen

Create the SpliceGraph of genes family given from a directory in hard code, the splicegraphs are in the spliceGraph directory
The paths should be put into the command line and not in hard code
Example : python spliceGraphMaker.py

'''

import sys
import requests
import os
from os import getcwd
from collections import Counter
import operator
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
from alternativeSplicing import *
from labelExons import *
from matrixScore import *
from clustering import *
from retrieveSeqClusters import *
from fasta_to_phylip import *
from microWithoutGenes import *
from compute_score_main import *
from fusionMatrix import *
from fuzzyCMeans import *
#from prepareTreeForConselTest import *

def build_arg_parser(path):
	parser = argparse.ArgumentParser(description="Ensembl gene tree parsor program parameters")   
	parser.add_argument('-s', '--save_path', default =  path + '/spliceGraph/')
	parser.add_argument('-ma', '--macroalignment', default = path + '/macroalignment/')
	parser.add_argument('-wst', '--weightStructure', default = 0.5)
	parser.add_argument('-wsq', '--weightSequence', default = 0.5)
	parser.add_argument('-i', '--inputFile', default = None)
	return parser

def one_spliceGraphMaker(save_path, root, filename):
	namefile = root +"/"+ filename
	exonslist = {}
	selectedexons = []
	nbtranscrit = 0
	numexon = 0	
	file = open(namefile)
	newname = filename.split("_")[0]
	namefile = newname + "_spliceGraphe.txt"
	completeName = os.path.join(save_path, namefile)         
	outputfile = open(completeName, "w")

	for line in file :
		#get the exon name 
		if line.startswith(">"):
			exonID = "exon_"
			exonID += str(numexon)
			numexon += 1
			exonslist[exonID] = []
		#when you read the genes and the transcripts
		else :

			#get the id 
			s = line.split(":")
			if len(s)>=2:
				i = s[0]
				l = s[1]

				#get the type
				t = i.split("0")
				if len(t)>=2:
					typ = t[0]
					val = t[1]
					if typ[-1] == "G":
						continue
					elif typ[-1] == "T":
						nbtranscrit += 1
						start = l.split("-")[0]
						end = l.split("-")[1]
						end = end.replace("\n", "")
						length = int(end) - int(start)
						position = [start, end]
						tab = [[start, end], length,i]
						exonslist[exonID].append(tab)	
						
	exonsID = exonslist.keys()

	for exon in sorted(exonsID):
		outputfile.write(">" + exon + "\n")
		positions = exonslist[exon]			
		locations = {}

		for position in positions:
			start = position[0][0]
			stop = position[0][1]
			length = position[1]
			transcritname = position[2]
			locationKey = start + "-" + stop
			if locationKey in locations.keys():
				locations[locationKey][0] += 1
			else:
				locations[locationKey] = [1,length, transcritname]

		#foreach exon get the one that is the most frequent in transcrits
		most_common = Counter(locations).most_common(nbtranscrit)
		
		#cannot loop on a 1 size array
		if len(most_common)-1 == 1:
			if most_common[0][0] == '0-0': #check if the exon is absent
				#check if the second exon is absent
				if most_common[1][0] != '0-0': 
					#if not then the higher value is the second exon
					outputfile.write(most_common[1][1][2] + ":" + str(most_common[1][0] + " \n")) 
			#if the first exon is present and if the frequency of first one is higher than the second than we take the first one
			elif most_common[0][1][0] > most_common[1][1][0]:
				outputfile.write(most_common[0][1][2] + ":" + str(most_common[0][0] + " \n"))
			#otherwise take the higher distance
			else: 
				l = (most_common[0][1][0],most_common[1][1][0])
				m = max(l)
				if l == most_common[0][1][0]:
					outputfile.write(most_common[0][1][2] + ":" + str(most_common[0][0]) + " \n")
				else: 
					outputfile.write(most_common[1][1][2] + ":" + str(most_common[1][0]) + " \n")
		
		#for array of bigger size
		if len(most_common)-1 >= 2:
			for k in range(len(most_common)-1):
				if most_common[k][0] == '0-0': 
					continue
				else:
					if most_common[k][1][0] > most_common[k+1][1][0]:
						outputfile.write(most_common[k][1][2] + ":" + str(most_common[k][0] + " \n"))
					else: #if fqz is equal
						if most_common[k][1][1] == most_common[k+1][1][1]: #if distance is equal
							outputfile.write(most_common[k][1][2] + ":" + str(most_common[k][0]) + " \n")
							#outputfile.write(most_common[k+1][1][2] + ":" + str(most_common[k+1][0]) + " \n")
						else: #distance is different, get max dist
							l = (most_common[k][1][1],most_common[k+1][1][1])
							m = max(l)
							if l == most_common[k][1][1]:
								outputfile.write(most_common[k][1][2] + ":" + str(most_common[k][0]) + " \n")
							else: 
								outputfile.write(most_common[k+1][1][2] + ":" + str(most_common[k+1][0]) + " \n")
					break					
	outputfile.close()

	return completeName, save_path, root, filename

def countNumberOfCDS(save_path, root,filename, macroalignment):	
	save_path = save_path.replace("spliceGraph", "initialSource")
	id = filename.split("_")[0]
	sourceToTarget = save_path + id + "_initialsource2target.txt"
	file = open(sourceToTarget, "r")
	lines = file.readlines()
	nb_cds = len(lines)

	if nb_cds <= 2:
		return False
	else:
		return True


def main_spliceGraphMaker(save_path, macroalignment,weightStructure, weightSequence):		
	for root, dirs, files in os.walk(macroalignment):  			
		ids = range(len(files))
		
		random.shuffle(ids)
		for i in ids:
			filename = files[i]		
		#for filename in files:	
			tmp = filename.split("_")
			prefixe= tmp[0]				
			if filename.count("~")>0:
				pass
				#break
			else:						
				#spliceGraphFile, save_path, root, filename = one_spliceGraphMaker(save_path, root, filename)							
				#executeAllScript(spliceGraphFile, save_path, root, filename)
				executeAllScript(save_path, root, filename, weightStructure, weightSequence)			
			

def spliceGraphMaker_aux(inputFilename):
	path =  sys.path[0]      
	save_path  = path + "/spliceGraph/"
	tmp = inputFilename.split("/")	
	filename = tmp[-1]
	tmp = tmp[:-1]
	root = "/".join(tmp)	
	spliceGraphFile, save_path, root, filename = one_spliceGraphMaker(save_path, root, filename)							
	executeAllScript(spliceGraphFile, save_path, root, filename)

def rewriteAln(input_file):
    records = SeqIO.parse(input_file, 'fasta')
    records = list(records) # make a copy, otherwise our generator
                            # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '-')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do alignment
    output_file = '{}.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    alignment = AlignIO.read(output_file, "fasta")
    
def executeAllScript(save_path, root, filename, weightStructure, weightSequence):
	try:
		minSimilarity = 0.5
		tmp = filename.split("_")
		prefixe= tmp[0]
		if countNumberOfCDS(save_path, root,filename, macroalignment) == False : #or os.path.exists("clusters_0_100/fuzzyCMeans/" + prefixe + "_cluster_1_NT.fasta"):
			return
		
		macroalignmentFile = root +filename
		microalignmentFile = path + "/microalignment/" + prefixe + "_microalignment.fasta"
		initialsource2target = path + "/initialSource/" + prefixe + "_initialsource2target.txt"

		rewriteAln(microalignmentFile)
		
		print("Compute labeledExons")
		labeledExons = main_labelExons_aux(macroalignmentFile, microalignmentFile)		


		print(labeledExons)
		print("remove Genes from microalignment")
		removedGeneInMicroALignment = microWithoutGenes(microalignmentFile)
		print (removedGeneInMicroALignment)
		print("Compute matrixScore based only on splice structure")
		matrixMicroalignment = main_matrixScore_aux(labeledExons, initialsource2target)	
		
		print("Compute matrixScore based on multiple alignment (Microalignment)")
		scoreMatrixSeq = main_compute_score(removedGeneInMicroALignment)

		print("Fusion of scoreMatrixSeq and matrixMicroalignment")
		for e in [[0.5,0.5]]:
		    weightStructure = e[0]
		    weightSequence = e[1]
		    sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))

		    if os.path.exists("/clusters_"+sub_path+"/fuzzyCMeans/" +  prefixe + "_cluster_1_NT.fasta") == True:
		        pass
		    else:
		        fusionSeqMat_MicroAlignMat, names, x_coords, u, mappingGeneToTranscript = main_fusion_matrix(scoreMatrixSeq, matrixMicroalignment, weightStructure, weightSequence)  

		        print("Compute clusters")
		        clustersFilenames = main_fuzzy_cMeans(names, x_coords, mappingGeneToTranscript, microalignmentFile, sub_path)


	except Exception as e:
		print(e)

		

if __name__ == "__main__":    

	path =  sys.path[0]
	parser = build_arg_parser(path)
	arg = parser.parse_args()        
	save_path  = arg.save_path 
	macroalignment = arg.macroalignment
	inputFile = arg.inputFile
	weightStructure = float(arg.weightStructure)
	weightSequence = float(arg.weightSequence)
	
	if inputFile == None:
		main_spliceGraphMaker(save_path, macroalignment, weightStructure, weightSequence) 		
	else:
		tmp = inputFile.split("/")		
		filename = tmp[-1]
		tmp = tmp[:-1]		
		root = "/".join(tmp)	
		tmp_savepath = save_path.split("/")
		tmp_savepath =tmp_savepath[-1]
		tmp_savepath = "/".join(tmp_savepath)
		executeAllScript(save_path, tmp_savepath+"macroalignment/", filename, weightStructure, weightSequence)

