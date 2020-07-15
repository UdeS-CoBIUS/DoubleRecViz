#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``utils.py`` **module description**:


This module defines functions required at all steps of the comparison.

.. moduleauthor:: Safa Jammali and Jean-David Aguilar
Université de Sherbrooke Canada
Laboratoty CoBiUS

2017-2018

"""

import os
import csv
from ete3 import *
from Bio import pairwise2


MAX_CDS_SEPARATION_FOR_CONCAT = 100
MAX_SEPARATION_DIFF_FOR_CONCAT = 60

def launch_tblastx(cdsfile, genefile,evalue, block1,block2,cdsid, geneid):
	"""
	This function 
	Parameters
	----------

	cdsfile:
	genefile:
	evalue:
	block1:
	block2:
	cdsid:
	geneid:

	Returns
	-------
	blastoutput:
	"""
	block1_qs,block1_qe,block1_ss,block1_se=block1
	block2_qs,block2_qe,block2_ss,block2_se=block2

	blastoutput = str(os.getcwd()+'/src/results/blast_results/CDS_'+ cdsid+'_'+str(block1_qe)+'-'+str(block2_qs)+'_VS_Gene_'+ geneid+'_'+str(block1_se)+'-'+str(block2_ss)+'.tblastx')
    
	command = "tblastx -query " + cdsfile + " -subject " + genefile +  " -evalue " + str(evalue) + " -outfmt "+ str(7) +  " -strand plus " +  " -out " + str(blastoutput)
	os.system(command)

	return blastoutput



def readBlastOut(blastOutput, blockQueryEnd, blockSubjectEnd,evalue):
	"""
	This function reads blast output file and extract the hits 

	Parameters
	----------
	bastOutput:
	blockQueryEnd:
	blockSubjectEnd

	Returns
	-------
	listHits:
	"""
	listHits=[]
	with open(blastOutput, 'rb') as csvfile:
		spamreader = csv.reader(csvfile,delimiter = ' ', quotechar = ' ')

		for row in spamreader:


			if row[0][0] is not '#':
				# construct a hit

				identity= float(row[0].split('\t')[2])
				beginHitCDS= int(row[0].split('\t')[6])+int(blockQueryEnd) - 1
				endHitCDS = int(row[0].split('\t')[7]) + int(blockQueryEnd)
				beginHitgene = int(row[0].split('\t')[8])+ int(blockSubjectEnd) - 1
				endHitgene = int(row[0].split('\t')[9]) + int(blockSubjectEnd)
				evalueHit = float(row[0].split('\t')[10])
				if (beginHitCDS < endHitCDS and beginHitgene < endHitgene and evalueHit <= evalue):
					listHits.append([beginHitCDS,endHitCDS, beginHitgene, endHitgene, evalueHit, identity])
	return listHits

	

def cover_percentage(blocklist,cdslength):
    """
    This function returns the percentage of coverage of the cds 
    by the blocklist.

    Parameters
    ----------

    blocklist:
    cdslength:

    Returns
    -------
    percentcover:
    """
    percentcover = 0.0
    
    coveredlength = 0

    for i in range(len(blocklist)):
        block = blocklist[i] 
        coveredlength += block[1]-block[0]
        if(i > 0):
            prev_block =  blocklist[i-1]
            if(block[2]-prev_block[3] == 0):
                coveredlength += block[0]-prev_block[1]
        
    percentcover= 100.0 * coveredlength / cdslength
    return percentcover

def order(blocklist):
    """
    This function orders blocks in the blocklist by increasing
    order of query start location.
    Parameters
    ----------

    blocklist:
    

    Returns
    -------
    blocklist:
   
    """

    for i in range(len(blocklist)):
        min = i
        for j in range(i+1,len(blocklist)):

            if(blocklist[j][0] < blocklist[min][0]):
                min = j
        if(min != i):
            tmp = blocklist[min]
            blocklist[min] = blocklist[i]
            blocklist[i] = tmp
    return blocklist


    


def computeAlignmentPercentIdentity(sequence1, sequence2):
    """
    This function aim to return the percentage of identity of 2 sequences
    sequence1: string
        sequence 1
    sequence2: string
        sequence 2

     Returns
     -------
     aln_identity: float
        % of identity between the 2 sequences
    """

    aln_identity = 0.0
    match = 0
    length = len(sequence1)

    for i in range(length):

        if (i < len(sequence2) and sequence1[i] == sequence2[i]):
            match += 1
    if length == 0:
        length = 1
    aln_identity = 100.0 * match / length
    return aln_identity
