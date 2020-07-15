#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``compute_score_main.py`` **module description**:

This module allows to compute the pairwise alignment scores of a multiple 
alignement of coding sequences given in a multifasta file and write the result 
score matrix in an output file. .

.. moduleauthor:: Aida Ouangraoua, Safa Jammali

April 2016

"""

import argparse
from Bio import SeqIO
import os
import time
import sys
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import numpy as np
import pandas as pd

LENGTH_ALN_LINE = 50
INFTY = 10000

def build_arg_parser(path):
    parser = argparse.ArgumentParser(description='Compute scores')
    parser.add_argument('-d', '--dataalignment')
    parser.add_argument('-o', '--outfile', default =  path + '/clusters/')
    return parser



def main():
    path =  sys.path[0]

    parser = build_arg_parser(path)

    arg = parser.parse_args()
    
    microalignmentFile = arg.dataalignment        

    main_compute_score(microalignmentFile)

def main_compute_score(microalignmentFile): 
    #rewriteAln(microalignmentFile)  
    tmp = microalignmentFile.split("/")   
    save_path = ""   
    for i in range(len(tmp) -2):
        save_path += tmp[i] + "/"
    """    
    filename = tmp[-1]
    tmp = tmp[:-1]      
    root = "/".join(tmp)
    path =  sys.path[0]

    parser = build_arg_parser(path)
    arg = parser.parse_args()   

    path =  sys.path[0] 

    save_path = path + "/similaritySeqScores/"
    """
    save_path +=  "similaritySeqScores/"
    
    microalignmentFilename = microalignmentFile.split("/")[-1]
    completeName = save_path + microalignmentFilename.split(".")[0] + ".matrix.csv"    

    aln = AlignIO.read(open(microalignmentFile), 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)

    names = dm.names
    n = len(names)

    distanceMatrix = dm.matrix

    matriceScore = np.eye(n, n)

    for seqId in names:
        distanceToOthers = dm[seqId]
        for j in range(n):
            i = names.index(seqId)
            matriceScore[i][j] = 1 - distanceToOthers[j]
            matriceScore[j][i] = 1 - distanceToOthers[j]

    

    df = pd.DataFrame(matriceScore, index=names, columns=names)

    with open(completeName, 'w') as outputfile:
        df.to_csv(outputfile)   


    outputfile.close()

    return completeName

if __name__ == '__main__':
    main()
