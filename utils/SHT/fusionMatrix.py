#!/usr/bin/python

'''
Authors : Esaie Kuitche and Marie Degen
Example : python -W ignore parseMatrixSeq.py similaritySeqScores/FAM86.matrix.txt similarityScores/FAM86_score.csv

'''

import sys
import argparse
import os
import ntpath
import numpy as np
import pandas as pd
from scipy.linalg import sqrtm
# generate related variables
from numpy import mean
from numpy import std
from numpy.random import randn
from numpy.random import seed
from matplotlib import pyplot
from numpy import cov
from scipy.stats import pearsonr
from skbio.stats.distance import mantel
import matplotlib.pyplot as plt
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


def normalizeScoreMatrix(scoreMatrix):
    scoreMatrix = scoreMatrix - scoreMatrix.min()
    scoreMatrix = scoreMatrix/float(scoreMatrix.max())
    return scoreMatrix;

def generateZero(number):
    """
    This function generate number zero
    :param number: the number of zero
    :return: number zero
    """
    stringOfZero = ""
    for i in range(number):
        stringOfZero = stringOfZero + "0\t"

    return stringOfZero

def provider(fileName):
    """
    This function take as entry a file containing a matrix, validate his format, correct his value because it's a symmetric matrix. After that it construct a numpy matrix that is nXn and return this numpy matrix
    :param file: the file name
    :return: numpy matrix in the form, and the columns and rows name
    """
    f = open(fileName, 'r')
    lines = f.readlines()
    header = ""
    stringOfMatrix = ""
    listeZero = ""
    i = 0
    for line in lines:
        line = line.strip()
        if (not (line.startswith("#")) and line != ""):
            if (i == 0):
                header = line.split("\t")
                i = i + 1
            else:
                listeZero = generateZero(i-1)
                i = i + 1
                line = listeZero + line + ";"
                stringOfMatrix = stringOfMatrix + line

    stringOfMatrix = stringOfMatrix #+ generateZero(i) + ";"

    #stringOfMatrix = stringOfMatrix[:-1]

    #print(stringOfMatrix)
    #xit()
    matrix = np.matrix(stringOfMatrix[:-1])
    return (matrix, header)
   
'''
Retrieve the value from the sequence matrix
'''
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

def main_normalize(matrixSequenceFile):
    matrixSequence, header = provider(matrixSequenceFile)
    matrixSequence = normalizeScoreMatrix(matrixSequence)
    return matrixSequence

def writeDistanceMartixInFil(distanceMatrix, header, matrixSequenceFileToSave):
    df = pd.DataFrame(distanceMatrix, index=header, columns=header)
    df.to_csv(matrixSequenceFileToSave, index=True, header=True, sep='\t')


def SVD(fusionMatrix, names):
    u, s, v = np.linalg.svd(fusionMatrix)    

    new_s = []
    for i in range(len(s)):
        if s[i]<1:
            break 
        else:
            new_s.append(s[i])

    n = len(new_s)

    new_u = u[:,:n]

    x_coords = {}

    #aa = np.dot(u,sqrtm(np.diag(s)))
    new_aa = np.dot(new_u,sqrtm(np.diag(new_s)))

    genes = {}

    transcriptName = []
    
    for name in names:
        tmp = name.split("_")
                
        if len(tmp)==2:
            transcript = tmp[0]
            gene = tmp[1]
        else:
            gene = tmp[-1]
            del tmp[-1]
            transcript = "_".join(tmp)       


        if gene in genes.keys():
            genes[gene].append(transcript)
        else:
            genes[gene] = [transcript]

        x_coords[transcript] = new_aa[names.index(name)]
        #x_coords[transcript] = u[names.index(name)]
        transcriptName.append(transcript)

        
    #print(new_u)
    #print(new_s)
    #print(transcriptName)
    #exit()
    return transcriptName, x_coords, new_aa, genes


def main_fusion_matrix(matrixSequenceFile, matrixStructureFile, weightStructure, weightSequence):
    path =  sys.path[0] 
    save_path = path + "/scores/"


    matrixStructure = pd.read_csv(matrixStructureFile)
    size = matrixStructure.shape[0]
    matrixStructureValues = matrixStructure.values

    #matrixSequence, header = provider(matrixSequenceFile)    
    
    matrixSequenceTmp = pd.read_csv(matrixSequenceFile)
    matrixSequenceTmp = matrixSequenceTmp.values

    header = retrieveName(matrixSequenceTmp)
    #print(header,size, matrixStructureFile,matrixSequenceFile)
    #exit()
    matrixSequence = retrieveValues(matrixSequenceTmp, size)     

    #exit()
    #matrixSequence = normalizeScoreMatrix(matrixSequence)

    matrixSequenceFileToSave = matrixSequenceFile.split(".")[0] + ".csv"

    writeDistanceMartixInFil(matrixSequence, header, matrixSequenceFileToSave)

    #exit()
    names = retrieveName(matrixStructureValues)
    values = retrieveValues(matrixStructureValues, size)  

    fusionMatrix = np.eye(size, size)
    data1 = []
    data2 = []
    for i in range(len(names)):
    	for j in range(len(names)):
            v = getValueSeqMatrix(i,j,matrixSequence, names, header)
            #data1.append(values[i][j])
            #data2.append(v)
            #print(i, j, v, values[i][j], float("{0:.2f}".format((values[i][j] + v)/2)))
            #print(i, j, values[i][j], v, values[i][j]-v )
            fusionMatrix[i][j] = fusionMatrix[j][i] = float("{0:.2f}".format( weightStructure*values[i][j] + weightSequence*v )) +0.000000000000000001
            #fusionMatrix[i][j] = fusionMatrix[j][i] = float("{0:.5f}".format(values[i][j]))    
    #covariance = cov(data1, data2)
    #print(covariance)    
    #corr, _ = pearsonr(data1, data2)
    #print('Pearsons correlation: %.3f' % corr)    
    #pyplot.scatter(data1, data2)
    #pyplot.show()
    namefile = matrixStructureFile
    namefile = namefile.split("/")[-1]
    namefile = namefile.split(".")[0]
    namefile = namefile + "_Final.csv"
    completeName = os.path.join(save_path, namefile)         

    df = pd.DataFrame(fusionMatrix, index=names, columns=names)

    with open(completeName, 'w') as outputfile:
    	df.to_csv(outputfile)

    transcriptName, x_coords, aa, genes = SVD(fusionMatrix, names)

    return completeName, transcriptName, x_coords, aa, genes

def correlation_matrix(matrixSequenceFile, matrixStructureFile):

    matrixStructure = pd.read_csv(matrixStructureFile)
    size = matrixStructure.shape[0]
    #print(size)
    if size <= 3:
        return [],[],[],[]

    matrixStructureValues = matrixStructure.values
    names = retrieveName(matrixStructureValues)
    #matrixSequence, header = provider(matrixSequenceFile)    
    
    matrixSequenceTmp = pd.read_csv(matrixSequenceFile)
    matrixSequenceTmp = matrixSequenceTmp.values

    header = retrieveName(matrixSequenceTmp)
    #print(header,size, matrixStructureFile,matrixSequenceFile)
    #exit()
    
    matrixSequenceValues = retrieveValues(matrixSequenceTmp, size)     
    transcriptName, x_coords, aa, genes = SVD(matrixSequenceValues, names)
    
    x_abs = []
    y_ord = []
    for x in transcriptName:
        if len(x_coords[x])>1:
            x_abs.append(x_coords[x][0])
            y_ord.append(x_coords[x][1])
            #plt.scatter(x_abs, y_ord, s=50, cmap='gray')
            #plt.show()

    matrixStructureValues = retrieveValues(matrixStructureValues, size) 
    transcriptName2, x_coords2, aa2, genes2 =  SVD(matrixStructureValues, names)
    x_abs2 = []
    y_ord2 = []
    for x in transcriptName2:
        if len(x_coords2[x])>1:
            x_abs2.append(x_coords2[x][0])
            y_ord2.append(x_coords2[x][1])
            #plt.scatter(x_abs, y_ord, s=50, cmap='gray')
            #plt.show()    
    #print(matrixSequenceValues)
    #print(matrixStructureValues)   
    """
    one_matrix = np.full((size, size), 1)
    print(matrixSequenceValues - matrixStructureValues)
    matrixSequenceValues =  (one_matrix - matrixSequenceValues)
    matrixStructureValues = (one_matrix - matrixStructureValues)
    
 
    correlation = mantel(matrixSequenceValues, matrixStructureValues, method='pearson')    
    """
    return x_abs, y_ord, x_abs2, y_ord2

def correlation_matrix2(matrixSequenceFile, matrixStructureFile):
    ratio = []
    matrixStructure = pd.read_csv(matrixStructureFile)
    size = matrixStructure.shape[0]


    matrixStructureValues = matrixStructure.values
    names = retrieveName(matrixStructureValues)
    #matrixSequence, header = provider(matrixSequenceFile)    
    
    matrixSequenceTmp = pd.read_csv(matrixSequenceFile)
    matrixSequenceTmp = matrixSequenceTmp.values

    
    matrixSequenceValues = retrieveValues(matrixSequenceTmp, size)     
 
    matrixStructureValues = retrieveValues(matrixStructureValues, size) 
 


    for i in range(len(names)):
        for j in range(len(names)):
            v1 = matrixSequenceValues[i,j]
            v2 = matrixStructureValues[i,j]
            if (v1!=0 and v2!=0):
                ratio.append(v1/v2)


    return ratio


def correlation_matrix3(matrixSequenceFile, matrixStructureFile):
    ratio = []
    matrixStructure = pd.read_csv(matrixStructureFile)
    size = matrixStructure.shape[0]
    if size < 10:
        return [],[]

    matrixStructureValues = matrixStructure.values
    names = retrieveName(matrixStructureValues)
    #matrixSequence, header = provider(matrixSequenceFile)    
    
    matrixSequenceTmp = pd.read_csv(matrixSequenceFile)
    matrixSequenceTmp = matrixSequenceTmp.values

    
    matrixSequenceValues = retrieveValues(matrixSequenceTmp, size)     
 
    matrixStructureValues = retrieveValues(matrixStructureValues, size) 
 
    X = []
    Y = []

    for i in range(len(names)):
        for j in range(i):
            #print(i,j)
            v1 = matrixSequenceValues[i,j]
            v2 = matrixStructureValues[i,j]
            if v1-v2>0.4 or v2-v1>0.4:
                pass
            else:
                X.append(v1)
                Y.append(v2)




    return X, Y
