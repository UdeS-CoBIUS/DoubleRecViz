	#!/usr/bin/python

'''
Author : Marie Degen


Cluster the transcrits from a genes family
Example : python -W ignore clustering.py similarityScores/FAM86_score.csv 0.5

'''

import sys
import argparse
import os
import numpy as np
import pandas as pd
import pickle
import json



def build_arg_parser(path):
	parser = argparse.ArgumentParser(description="Label exons")   
	parser.add_argument('-s', '--save_path', default =  path + '/clusters/')
	parser.add_argument('-i', '--inputFile',  required=True)
	parser.add_argument('-m', '--minSimilarity',  default = 0.5)
	return parser



'''
Return true if the transcrits are from the same gene, false otherwise
'''
def checkSameGene(transcrit1, transcrit2):	
	if len(transcrit1.split("-")) < 2  and len(transcrit2.split("-")) < 2:
		return True if transcrit1.split("_")[1] == transcrit2.split("_")[1] else False
	elif len(transcrit1.split("-")) > 2  and len(transcrit2.split("-")) < 2:
		for gene in transcrit1.split("-"):
			if gene.split("_")[1] == transcrit2.split("_")[1]:
				return True
			else: 
				continue 
	elif len(transcrit1.split("-")) < 2  and len(transcrit2.split("-")) > 2:
		for gene in transcrit2.split("-"):
			if gene.split("_")[1] == transcrit1.split("_")[1]:
				return True
			else: 
				continue 
	elif len(transcrit1.split("-")) > 2  and len(transcrit2.split("-")) > 2:
		for gene1, gene2 in zip(transcrit2.split("-"), transcrit1.split("-")):
			if gene1.split("_")[1] == gene2.split("_")[1]:
				return True
			else: 
				continue 

'''
Return true if the transcrits are the same
'''
def checkSameTranscrit(transcrit1, transcrit2):
	return True if transcrit1.split("_")[0] == transcrit2.split("_")[0] else False

'''
Retrieve the names of the header
'''
def retrieveName(matrix):
	names = []
	for val in range(len(matrix)): 
		names.append(matrix[val][0].replace("\n",""))
	return names

'''
Reconstruct the matrix from the file, return a matrix with only the values and a dataframe with the header
'''
def retrieveValues(matrixValues, size):	
	values = np.eye(size, size+1)

	for i, v in enumerate(matrixValues): 
		for j, o in enumerate(v): 
			if type(o) == str  : 
				continue 
			else:
				values[i][j] = o
	
	values = np.delete(values,0,1)
	return values
	
'''
Return the maximum value of similarity found in the matrix, the two transcrits corresponding
'''
def returnMaxValueSimilarity(matrix, names, bannedGroup):	
	maxValue = 0
	tr1 = ""
	tr2 = ""

	for i in range(len(matrix)):
		for j in range(len(matrix)):		
			if names[i] == names[j]:
				continue
			else: 				
				if checkSameGene(names[i],names[j]):
					continue
				elif names[i].replace("\n","") in bannedGroup or names[j].replace("\n","") in bannedGroup: 
					continue
				else:					
					if matrix[i][j] > maxValue: 
						tr1 = names[i]
						tr2 = names[j]
						maxValue = matrix[i][j]
	return maxValue, tr1, tr2

'''
Return the matrix to join to the bigger matrix, the matrix to join has the average calculated in it 
'''
def matrixDistCalculate(matrix, transcrit1, transcrit2, bannedGroup, names):
	matrixToJoin = []
	for k in range(len(matrix)):
		nbvalue = 0
		moyenne = 0 
		summoyenne =0 
		for j in range(len(matrix)):
			if checkSameGene(names[k],transcrit2) or checkSameGene(names[k],transcrit1) :
				if len(transcrit2.split("-")) > 2 or len(transcrit1.split("-")) > 2:
					if checkSameGene(names[k],transcrit2) or checkSameGene(names[k],transcrit2):
						continue
					elif checkSameGene(names[k],transcrit1) or checkSameGene(names[k],transcrit1):
						continue
				else:
					continue
			elif names[k] in bannedGroup:
				continue
			else:
				summoyenne += matrix[k][j]
				nbvalue += 1
			try: 
				moyenne = summoyenne / nbvalue
			except:
				moyenne = 0.0
		matrixToJoin.append(float("{0:.2f}".format(moyenne)))
	matrixToJoin.append(1.0)
	return matrixToJoin


def main_cluster(save_path, inputFile, minSimilarity):
	size = 0
	matrixValues = []
	cluster = []
	bannedGroup = []
	initGroup = []

	matrixScore = pd.read_csv(inputFile)
	size = matrixScore.shape[0]
	matrixValues = matrixScore.values

	names = retrieveName(matrixValues)
	values = retrieveValues(matrixValues, size)
	initGroup = names[:]
	namefiletmp = inputFile.split("/")[-1]
	namefile = namefiletmp.split(".")[0]
	#familyName = namefiletmp.split("_")[1]
	#namefile = namefile + "_" + familyName + "_clusters.out"
	namefile = namefile + "_clusters.out"
	completeName = os.path.join(save_path, namefile)         
	outputfile = open(completeName, "w")

	while True:
		maximum, transcrit1, transcrit2 = returnMaxValueSimilarity(values, names, bannedGroup)	
		if maximum <= float(minSimilarity):
			break
		else:
			#delete the transrit from the finalgroup
			initGroup.remove(transcrit1.replace("\n",""))
			initGroup.remove(transcrit2.replace("\n",""))
			
			name_cluster = "(" + transcrit1.replace("\n","") + "-" + transcrit2.replace("\n","") + ")"
			
			cluster.append([transcrit1.replace("\n",""), transcrit2.replace("\n",""), maximum])
			bannedGroup.append(transcrit1.replace("\n",""))
			bannedGroup.append(transcrit2.replace("\n",""))
			names.append(name_cluster)
			initGroup.append(name_cluster)
			
			matrixToJoin = matrixDistCalculate(values, transcrit1, transcrit2, bannedGroup, names)
			
			#resize the matrix with the new values
			n = values.shape[0]
			w = np.eye(n+1, n+1)
			w[:-1,:-1] = values
			w[-1:] = matrixToJoin
			w[:,-1] = matrixToJoin
			values = w
				
	outputfile.write(str(initGroup))
	return completeName
	#pickle.dump(initGroup, outputfile)
	#json.dump(initGroup, outputfile)			


def main_cluster_aux(inputFile, minSimilarity):
	path =  sys.path[0]      
	save_path  = path + "/clusters/"
	return main_cluster(save_path, inputFile, minSimilarity)


if __name__ == "__main__":    
	path =  sys.path[0]
	parser = build_arg_parser(path)
	arg = parser.parse_args()        
	save_path  = arg.save_path 
	inputFile = arg.inputFile
	minSimilarity = arg.minSimilarity
	main_cluster(save_path, inputFile, minSimilarity)
