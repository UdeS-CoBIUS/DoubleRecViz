#!/usr/bin/python

'''
Author : Marie Degen


Calculate the difference between the matrix generated from the two methods 
Example : python -W ignore diffmatrix.py similarityScores/FAM86_score.csv similarityScores/FAM86_score1.csv

'''

import sys
import argparse
import os
import numpy as np
import pandas as pd
import pickle
import json

#Save Path 
save_path = '/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph'

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
args = parser.parse_args()

mat1 = pd.read_csv(sys.argv[1])
size = mat1.shape[0]
mat1Values = mat1.values

mat2 = pd.read_csv(sys.argv[2])
mat2Values = mat2.values

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
def retrieveValues(matrix):
	values = np.eye(size, size+1)
	for i, v in enumerate(matrix): 
		for j, o in enumerate(v): 
			if type(o) == str  : 
				continue 
			else:
				values[i][j] = o
	
	values = np.delete(values,0,1)
	return values

names = retrieveName(mat1Values)
values = retrieveValues(mat1Values)  

names2 = retrieveName(mat2Values)
values2 = retrieveValues(mat2Values)  


diffMatrix = np.eye(values.shape[0], values.shape[0])
for i in range(len(names)):
	for j in range(len(names)):
		diffMatrix[i][j] = diffMatrix[j][i] = float("{0:.2f}".format(values[i][j] - values2[i][j]))
		
df = pd.DataFrame(diffMatrix, index=names, columns=names)

with open("diffmatrix.csv", 'wb') as outputfile:
	df.to_csv(outputfile)

print diffMatrix
		
		
		
		
		
		
		
		
		
