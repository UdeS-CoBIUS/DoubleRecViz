 #!/usr/bin/python
import sys
import requests
import os
from os import getcwd
import argparse
from ete3 import Tree
import glob
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
import pandas as pd
import random, string
from Bio.Align.Applications import MafftCommandline
#from StringIO import StringIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from copy import deepcopy, copy
import traceback
import numpy as np
from scipy.stats.stats import pearsonr  
import math
import matplotlib.pyplot as plt
import pandas as pd

def main_corrected_ensemblTree(all=0):	
	SHT_similaritry = [[],[],[],[],[],[],[]]
	LT_similaritry = [[],[],[],[],[],[],[]]
	Ration_SHT_LT_similaritry = [[],[],[],[],[],[],[]]
	for x in glob.glob("clusters_100_0/fuzzyCMeans/*_cluster_1_root.nhx"): 		
	    try:
		    x2 =x
		    file = open(x, "r")
		    x = str(x)
		    x= x.split("/")
		    x = x[-1]		
		    x = x.split("_")[0]
		    #print(x)
		    mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
		    mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")

		    mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")

		    trees_build_by_orthoGroup_nw = x2
		    mapping_genes = {}
		    lines = mapping_genes_file.readlines()	
		    #print(lines)
		    for line in lines:
			    line = line.replace("\n", "")
			    line = line.split("\t")
			    mapping_genes[line[0]] = [line[1], line[2], line[3]]
			

		    file = open(trees_build_by_orthoGroup_nw, "r")
		    tree_str = file.read()
		    tree_str = tree_str.replace("\n", "")
		    tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
		    tree = Tree(tree_str)
		    leavesNames = []
		    SHT_CDS_dict = {}
		    SHT_CDS = []

		    for node in tree.traverse("postorder"):
			    if node.is_leaf():
				    name = node.name
				    SHT_CDS_dict[mapping_genes[name][1]] = mapping_genes[name][0]
				    name =  mapping_genes[name][0]
				    SHT_CDS.append(name)

		    #print("isoselSeq/" + x + "_root.nhx")
		    file = open("isoselSeq/" + x + "_root.nhx" , "r")
		    tree_str = file.read()
		    tree_str = tree_str.replace("\n", "")
		    tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
		    tree = Tree(tree_str)

		    leavesNames = []
		    LT_CDS = []
		    LT_CDS_dict = {}
		    genes_used_by_LT = []

		    #print(SHT_CDS)
		    #print(mapping_genes)
		    for node in tree.traverse("postorder"):
			    if node.is_leaf():
				    name = node.name
				    if name in mapping_genes.keys():

				        var = mapping_genes[name][1]
				        LT_CDS_dict[mapping_genes[name][1]] = mapping_genes[name][0]

				        name =  mapping_genes[name][0]

				        LT_CDS.append(name)

				        genes_used_by_LT.append(var)



		    #print(LT_CDS)
		    #print(SHT_CDS)
		    for k, v in SHT_CDS_dict.items():
			    if k not in genes_used_by_LT:
				    SHT_CDS.remove(v)

		    values = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1]
		    if all == 1 and len(LT_CDS)==len(SHT_CDS):
			    for alpha in values:
				    sim_LT, sim_SHT = mean_simalarity(x, LT_CDS, SHT_CDS, alpha)

				    SHT_similaritry[values.index(alpha)].append(sim_SHT)
				    LT_similaritry[values.index(alpha)].append(sim_LT)
				    Ration_SHT_LT_similaritry[values.index(alpha)].append(SHT_similaritry/LT_similaritry)
				    print(x, len(LT_CDS), SHT_similaritry/LT_similaritry)

		    elif all ==0 and len(LT_CDS)==len(SHT_CDS):
			    if (set(SHT_CDS) != set(LT_CDS)):
				    for alpha in values:
					    sim_LT, sim_SHT = mean_simalarity(x, LT_CDS, SHT_CDS, alpha)

					    SHT_similaritry[values.index(alpha)].append(sim_SHT)
					    LT_similaritry[values.index(alpha)].append(sim_LT)
					    Ration_SHT_LT_similaritry[values.index(alpha)].append(sim_SHT/sim_LT)
					    print(x, len(LT_CDS), sim_SHT/sim_LT)
	    except Exception as e:
	        #print(e)
	        pass
	#print(SHT_similaritry)
	#print(LT_similaritry)
	c = "black"
	plt.boxplot(Ration_SHT_LT_similaritry,
	                     vert=True,  # vertical box alignment
	                     #patch_artist=True,  # fill with color
	                      boxprops=dict(color=c),
                          capprops=dict(color=c),
                          whiskerprops=dict(color=c),
                          flierprops=dict(color=c, markeredgecolor=c),
                          medianprops=dict(color=c),
	                     labels=[r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks
	plt.xticks(fontsize=70, rotation=50)
	plt.yticks(fontsize=70)
	plt.ylabel('Ratio of mean similarity', fontsize=60)
	plt.legend(prop={'size': 60})	
	plt.show()

	
	
	mean_std_SHT = [[np.mean(x), np.std(x)] for x in SHT_similaritry]
	mean_std_LT = [[np.mean(x), np.std(x)] for x in LT_similaritry]
	mean_std_ratio = [[np.mean(x), np.std(x)] for x in Ration_SHT_LT_similaritry]
	print(mean_std_SHT)
	print(mean_std_LT)
	print(mean_std_ratio)

	mean_SHT = [item[0] for item in mean_std_SHT]
	std_SHT = [item[1] for item in mean_std_SHT]

	mean_LT = [item[0] for item in mean_std_LT]
	std_LT = [item[1] for item in mean_std_LT]	

	mean_ratio = [item[0] for item in mean_std_ratio]
	std_ratio = [item[1] for item in mean_std_ratio]	

	print(mean_SHT)
	print(std_SHT)
	print(mean_LT)
	print(std_LT)
	print(mean_ratio)	
	print(std_ratio)
def retrieveName(matrix):
	names = []
	for val in range(len(matrix)): 
		names.append(matrix[val][0].replace("\n",""))
	return names

'''
Reconstruct the matrix from the file, return a matrix with only the values and a dataframe with the header
'''
def retrieveValues(matrix, size):
	values = np.eye(size, size+1)
	for i, v in enumerate(matrix): 
		for j, o in enumerate(v): 
			if type(o) == str  : 
				continue 
			else:
				values[i][j] = o
	
	values = np.delete(values,0,1)
	return values
sht = ["ENSAHAT00000008308","ENSAOWT00000032439","ENSARWT00000019779","ENSDNVT00000007356","ENSNPET00000021093","ENSCPGT00000028508","ENSCPUT00000030256","ENSZALT00000005102"]
lt = ["ENSZALT00000005102","ENSNPET00000021018","ENSFALT00000004623","ENSDNVT00000007356","ENSCPUT00000030252","ENSCPGT00000028560","ENSARWT00000019779","ENSAOWT00000032439","ENSAHAT00000008308"]

def mean_simalarity(x, LT_CDS, SHT_CDS, alpha):
	struc_file = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/similarityScores/microalignment_"+x+"_score.csv"
	seq_file = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/similaritySeqScores/"+x+"_microalignment.csv"
	matrixStruct = pd.read_csv(struc_file)		
	sizestruc = matrixStruct.shape[0]
	matrixStructValues = matrixStruct.values
	namesStruct = retrieveName(matrixStructValues)
	#namesStruct = [x.split("_")[0] for x in namesStruct]
	matrixStructValues = retrieveValues(matrixStructValues, sizestruc)	 


	matrixSeq = pd.read_csv(seq_file, sep="\t")
	sizeseq = matrixSeq.shape[0]
	matrixSeqValues = matrixSeq.values
	namesSeq = retrieveName(matrixSeqValues)
	matrixSeqValues = retrieveValues(matrixSeqValues, sizeseq) 
	LT_sim = 0.0
	SHT_sim = 0.0
	cmpt1 = 0
	#print(LT_CDS)
	#print(SHT_CDS)
	for i in range(len(namesStruct)):
		for j in range(i):
			#print(namesStruct[i].split("_")[0], namesStruct[j])
			if (namesStruct[i].split("_")[0] in LT_CDS) and  (namesStruct[j].split("_")[0] in LT_CDS):
				cmpt1 +=1
				stuct_sim = float(matrixStructValues[i,j])
				seq_sim = float(getValueSeqMatrix(i,j,matrixSeqValues, namesStruct, namesSeq))
				LT_sim += stuct_sim*alpha + (1-alpha)*seq_sim
	cmpt2 = 0
	#print("______________________________--")
	for i in range(len(namesStruct)):
		for j in range(i):
			#print(namesStruct[i], namesStruct[j])
			if (namesStruct[i].split("_")[0] in SHT_CDS) and  (namesStruct[j].split("_")[0] in SHT_CDS):
				cmpt2 +=1
				stuct_sim = float(matrixStructValues[i,j])
				seq_sim = float(getValueSeqMatrix(i,j,matrixSeqValues, namesStruct, namesSeq))
				SHT_sim += stuct_sim*alpha + (1-alpha)*seq_sim

	if cmpt1 ==0 or cmpt2 == 0:
		#exit()
		return -1, -1
	else:
		return LT_sim/cmpt1, SHT_sim/cmpt2

def getValueSeqMatrix(position1, position2, matrixSequence, names, header):
	namesonly = []

	for nameandgene in names:
		tmp = nameandgene.split("_")
		tmp_name = ""
		
		if len(tmp)==2:
			namesonly.append(tmp[0])
		else:
			del tmp[-1]
			tmp_name = "_".join(tmp)		
			namesonly.append(tmp_name)


	for i, position in enumerate(header):
		if position == namesonly[position1]:
			pos1 = i

	for j,pos in enumerate(header): 
		if pos == namesonly[position2]:
			pos2 = j


	return matrixSequence[pos1, pos2]


def correlation_matrix():
	correlations = []
	correlations2 = []
	for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/similaritySeqScores/*_microalignment.csv"): 
		seq_file =x
		file = open(x, "r")
		x = str(x)
		x= x.split("/")
		x = x[-1]		
		x = x.split("_")[0]
		print(x)
		struc_file = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/similarityScores/microalignment_" + x + "_score.csv"
		if os.path.exists(struc_file):
			matrixStruct = pd.read_csv(struc_file)		
			sizestruc = matrixStruct.shape[0]
			matrixStructValues = matrixStruct.values
			namesStruct = retrieveName(matrixStructValues)
			#namesStruct = [x.split("_")[0] for x in namesStruct]
			matrixStructValues = retrieveValues(matrixStructValues, sizestruc)	 


			matrixSeq = pd.read_csv(seq_file, sep="\t")
			sizeseq = matrixSeq.shape[0]
			matrixSeqValues = matrixSeq.values
			namesSeq = retrieveName(matrixSeqValues)
			matrixSeqValues = retrieveValues(matrixSeqValues, sizeseq) 
			struc_values = []
			seq_values = []
			flag = False
			for e in namesStruct:
				if len(e.split("_")) > 2:
					flag = True
			if flag == False:				
				for i in range(len(namesStruct)):
					for j in range(i):
						stuct_sim = float(matrixStructValues[i,j])
						seq_sim = float(getValueSeqMatrix(i,j,matrixSeqValues, namesStruct, namesSeq))
						struc_values.append(stuct_sim)
						seq_values.append(seq_sim)
				correlation = pearsonr(struc_values, seq_values)
				
				correlation2 = np.corrcoef(struc_values,seq_values)

				if np.isnan(correlation[0]) or correlation[0]<0:
					#print(struc_values, seq_values)
					pass
				else:
					correlations.append(correlation[0])

				if np.isnan(correlation2[1][0]) or correlation2[1][0]<0.2:
					#print(struc_values, seq_values)
					pass
				else:
					correlations2.append(correlation2[1][0])

	
	#print(correlations)
	print(np.mean(correlations))
	print(np.median(correlations))
	
	#print(correlations2)
	print(np.mean(correlations2))
	print(np.median(correlations2))
	print(np.min(correlations2), np.max(correlations2))
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)
	plt.legend(prop={'size': 25})

	df = pd.DataFrame({'Correlation':correlations2})
	df.boxplot()
	plt.show()


def shareleaves():
	shares1 = [66.67,100,90,83.33,62.5,61.11,75,33.33,56.25,100,76.47,100,58.33,100,100,100,66.67,87.5,80,60,100,100,100,50,50,100,76.92,66.67,33.33,66.67,100,100,100,100,75,100,87.5,100,62.5,100,100,57.14,100,83.33,80,100,50,90.91,80,50,90,100,94.74,90,100,87.5,28.57,100,72.73,71.43,100,84.21,91.67,100,100,66.67,80,100,53.33,100,100,88.89,57.14,69.23,75,100,100,100,100,100,70,72.73,66.67,100,100,88.24,78.57,100,66.67,75,75,57.14,100,55.56,100,80,40,92,0,100,90,81.4437623762,17.0092932065]
	shares2 = [66.67,100,90,50,87.5,88.89,62.5,100,50,100,76.47,100,58.33,100,100,100,66.67,87.5,80,100,100,100,75,50,46.67,100,76.92,77.78,66.67,100,100,100,100,80,68.75,100,75,100,62.5,100,100,71.43,100,83.33,80,100,75,100,80,75,80,100,84.21,80,100,87.5,33.33,28.57,100,63.64,71.43,100,78.95,91.67,100,100,66.67,60,66.67,40,94.12,100,55.56,57.14,69.23,33.33,75,100,100,100,100,100,50,80,63.64,66.67,100,85.71,88.24,85.71,100,100,75,75,71.43,100,55.56,100,80,40,96,100,100,81.9277669903,16.0442831558]
	shares3 = [100,100,66.67,62.5,77.78,75,50,100,88.24,100,66.67,100,100,100,100,87.5,80,100,77.78,100,50,33.33,53.33,100,76.92,33.33,66.67,83.33,100,100,81.25,87.5,87.5,100,87.5,100,100,71.43,100,83.33,80,100,50,90.91,60,56.25,90,100,89.47,90,100,87.5,28.57,100,63.64,71.43,100,78.95,75,100,100,66.67,80,100,80,100,100,77.78,57.14,76.92,44.44,75,100,100,100,100,80,63.64,83.33,85.71,88.24,64.29,100,66.67,75,75,71.43,100,55.56,100,60,40,96,100,90,82.0221052632,15.3875900277]
	shares4 = [33.33,100,90,50,75,88.89,62.5,33.33,37.5,100,82.35,100,83.33,100,100,100,62.5,80,60,100,88.89,100,75,50,40,100,76.92,66.67,33.33,100,100,100,100,87.5,100,100,88.89,75,100,100,100,100,83.33,80,100,75,100,37.5,80,100,89.47,90,100,87.5,28.57,100,81.82,92.86,100,78.95,91.67,100,100,66.67,70,66.67,66.67,100,100,66.67,57.14,69.23,44.44,100,100,100,100,100,66.67,90,63.64,85.71,94.12,85.71,100,66.67,50,75,57.14,100,55.56,100,60,40,92,100,90,81.4155670103,17.3595281114]
	shares5 = [66.67,84.62,100,100,66.67,87.5,72.22,75,66.67,37.5,100,70.59,100,83.33,100,100,100,100,62.5,60,60,60,100,88.89,100,25,50,60,100,84.62,88.89,33.33,0,100,100,100,80,87.5,87.5,87.5,88.89,75,100,100,100,100,83.33,80,100,75,100,40,62.5,80,100,89.47,90,100,87.5,33.33,28.57,100,72.73,85.71,100,52.63,91.67,100,100,100,80,66.67,73.33,100,100,44.44,57.14,69.23,55.56,100,100,100,100,100,66.67,70,63.64,85.71,94.12,85.71,100,66.67,25,75,57.14,100,55.56,100,60,40,96,100,90,79.8099029126,18.2608954661]
	shares6 = [33.33,100,90,50,75,94.44,62.5,66.67,31.25,100,88.24,100,83.33,100,100,100,100,60,60,100,88.89,100,25,50,53.33,100,84.62,88.89,66.67,66.67,100,100,100,100,100,87.5,87.5,88.89,75,100,100,100,100,83.33,66.67,100,50,100,40,62.5,80,100,94.74,100,100,100,33.33,28.57,100,81.82,85.71,100,78.95,100,100,100,66.67,60,100,60,94.12,100,88.89,57.14,69.23,44.44,100,100,100,100,100,90,63.64,100,100,94.12,78.57,100,66.67,50,75,57.14,100,55.56,100,60,40,96,100,100,82.1053,18.373628]
	shares7 = [66.67,84.62,100,90,66.67,87.5,55.56,62.5,66.67,31.25,100,76.47,100,83.33,100,100,100,100,87.5,60,60,100,88.89,100,25,50,60,100,84.62,100,66.67,66.67,100,100,100,80,100,100,87.5,88.89,62.5,100,100,85.71,100,83.33,66.67,100,50,100,81.25,80,100,94.74,80,100,100,33.33,42.86,100,81.82,85.71,100,84.21,100,100,100,66.67,90,100,73.33,100,100,88.89,57.14,76.92,100,100,100,100,100,100,80,63.64,100,100,100,85.71,100,66.67,50,75,71.43,100,55.56,100,60,40,96,100,90,83.9214851485,15.7892657583]
	c = "black"
	plt.boxplot([shares1, shares2, shares3, shares4, shares5, shares6, shares7],
	                     vert=True,  # vertical box alignment
	                      boxprops=dict(facecolor=c, color=c),
                          capprops=dict(color=c),
                          whiskerprops=dict(color=c),
                          flierprops=dict(color=c, markeredgecolor=c),
                          medianprops=dict(color=c),
	                     labels=[r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks
	plt.xticks(fontsize=100, rotation=50)
	plt.yticks(fontsize=70)
	plt.ylabel('Ratio of mean similarity', fontsize=60)
	plt.legend(prop={'size': 60})
	plt.show()

#print(mean_simalarity("ENSGT00940000155862", lt, sht, 1))
#exit()
#correlation_matrix()		
main_corrected_ensemblTree()		
#shareleaves()
