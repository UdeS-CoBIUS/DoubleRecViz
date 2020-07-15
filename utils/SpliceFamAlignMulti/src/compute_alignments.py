#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compute_alignments.py`` **module description**:

This module is the module that launches the pairwise comparison for all pairs of source CDS and target gene

moduleauthor:: Safa Jammali
Université de Sherbrooke Canada

2017-2018

"""

import multiprocessing
from functools import partial
#from contextlib import contextmanager
from multiprocessing import Pool

import os
from utils import *
from step0_structure import *
from step1_local import *
from step2_extend import *
from  step3_global import *
from write import *

# Constant Declaration
STATUS_EXISTING_PROTEIN = 1
STATUS_PREDICTED_PROTEIN = 2
STATUS_COMPLETE_CDS = 3
STATUS_COMPLETE_CDS_DIFFERENT_STRUCTURE = 4
STATUS_PARTIAL_CDS = 5
EVALUE_STRUCTURE= 1.0/10000000

CANONICAL_DONOR = "GT"
CANONICAL_ACCEPTOR = "AG"
NC_DONOR1="GT"
NC_ACCEPTOR1="AG"
NC_DONOR2="GT"
NC_ACCEPTOR2="AG"

#########################
### COMPARISON ##########
#########################


def spliceAlignment(sourceData, targetData,  step, compareExon, choice, pairwiseMethod, outputPrefix):

    """
    This function launches a parallel pairwise alignment between all pairs of CDS and gene

    Parameters
    ----------
    sourceData: list
    	list that contains information about CDS : CDS ID, CDS sequence, gene ID, and exon list
    targetData: list
    	list that contains information about genes: gene id, and gene sequence
    step: string (1, 2, 3)
    	3: execute the DP step (Step3) of the method ;  2 : stop after the extension step ; 1: stop after local alignment
    choice: string
	splign ;  blast
    outputPrefix:string
	prefix for the output 

    Returns
    -------
    comparisonresults: list
        matrix of lists of alignment blocks (comparisonresults[i][j] : alignments of ith gene and jth CDS)
    """

    comparisonResults = []
    exonEndList = {}
    exonStartList = {}
    cdsExon = {}
    cds2GeneExon = {}
    geneExon = {}
    intronList = {}


    # put each cds in a fasta format file
    for cds in sourceData:
        cdsId, cdsSeq, null, null = cds
        cdsSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
        cdsFile = open(cdsSeqFile, "w")
        cdsFile.write(">" + cdsId + "\n")
        cdsFile.write(cdsSeq)
        cdsFile.close()

    # put each gene in a fasta format file
    for gene in targetData:
        geneId, geneSeq = gene
        geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId + '.fasta'
        geneFile = open(geneSeqFile, "w")
        geneFile.write(">" + geneId + "\n")
        geneFile.write(geneSeq)
        geneFile.close()


    # initialize tables for gene splicing structure
    for gene in targetData:
        geneId, geneSeq = gene
        exonStartList[geneId] = []
        exonEndList[geneId] = []
        geneExon[geneId] = []
        intronList[geneId] = []

        
    # extract splicing structure information for each gene based on its CDS
    sourceData2=[]
    for cds in sourceData:
	    blockList_ =[]
            cdsId, cdsSeq, cdsGeneId, blockList = cds
            cdsExon[cdsId] = []
            cdsSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
            geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + cdsGeneId  + '.fasta'
            # if no splicing sructure information was given in input
            if(len(blockList) == 0):
		if  choice == 'splign':
		        splignOutputName = launch_splign(cdsSeqFile, geneSeqFile, cdsId, cdsGeneId )
		        blockList_ = parse_splign_output(splignOutputName)
			blockList = blockList_
                       
		elif choice == 'blast':
                    blockList_, blockList = BlastLocal(cdsSeqFile, geneSeqFile,  cdsId, cdsGeneId, EVALUE_STRUCTURE)
		else:
			print 'choose method to extractexon strusture by splign or blast'
           		exit(-1)
            # if the splicing sructure information was given
            else:
                blockList_ = blockList
                        
	    sourceData2.append([cdsId, cdsSeq, cdsGeneId, blockList_])

            # fill tables for gene splicing structure
            intronList[cdsGeneId ] += extractIntronGene(cdsGeneId , cdsId, blockList)
            exonEndList[cdsGeneId ] += [block[3] for block in blockList]
            exonStartList[cdsGeneId ] += [block[2] for block in blockList]
            geneExon[cdsGeneId ] += [block[2:] for block in blockList]
            cdsExon[cdsId] = [block[:2] for block in blockList]
            checkExon(cdsId, cdsExon[cdsId])
            cds2GeneExon[cdsId] = [block[2:] for block in blockList]
   
    if len(sourceData2[0][3])!=0:
	    write_input_files(outputPrefix,sourceData2,targetData)
    
    # Order gene splicing structure information
    for gene in targetData:
        geneId, geneSeq = gene
        exonStartList[geneId].sort()
        exonEndList[geneId].sort()
        geneExon[geneId].sort()
        intronList[geneId].sort()


        # remove exon and intron redundancy
        for i in range(len(geneExon[geneId]) - 1, 0, -1):
            if (geneExon[geneId][i] == geneExon[geneId][i - 1]):
                geneExon[geneId].remove(geneExon[geneId][i])
        for i in range(len(intronList[geneId]) - 1, 0, -1):
            if (intronList[geneId][i] == intronList[geneId][i - 1]):
                intronList[geneId].remove(intronList[geneId][i])
        for i in range(len(exonStartList[geneId]) - 1, 0, -1):
            if (exonStartList[geneId][i] == exonStartList[geneId][i - 1]):
                exonStartList[geneId].remove(exonStartList[geneId][i])
        for i in range(len(exonEndList[geneId]) - 1, 0, -1):
            if (exonEndList[geneId][i] == exonEndList[geneId][i - 1]):
                exonEndList[geneId].remove(exonEndList[geneId][i])

    # lauch pairwise comparison for all pairs of gene and CDS
    gene_cds = []
    for gene in targetData:
        for cds in sourceData2:
            gene_cds.append([gene,cds])
    p = Pool(multiprocessing.cpu_count())
    #pool_results = p.map(partial(poolCompare, cdsExon = cdsExon, sourceData = sourceData, cds2GeneExon = cds2GeneExon, intronList = intronList, geneExon = geneExon,  force = force), gene_cds)

    pool_results = []
    for x in gene_cds:
        pool_results.append(poolCompare(x,cdsExon,  sourceData, cds2GeneExon,intronList,geneExon,step, compareExon,pairwiseMethod))
    
    results = {}
    # recover results of all pairs
    for item in pool_results:
        results[item[5] + item[6]] = [item[0], item[1], item[2], item[3], item[4]]
   
    # delete CDS fasta files created before
    for cds in sourceData:
        cdsId, cdsSeq, null, null = cds
        cdssSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
        os.remove(cdssSeqFile)
    
    # delete gene fasta files and put results in a list according to the order of gene and CDS
    for gene in targetData:
        geneId, geneSeq = gene
        geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId + '.fasta'
        os.remove(geneSeqFile)
        comparisonResults.append([])
        for cds in sourceData:
            cdsId, null, null, null = cds
            comparisonResults[-1].append(results[geneId + cdsId])
    
    # return matrix of lists of blocks of alignments
    return comparisonResults


def checkExon(cdsId, blockList):
    '''
    This function check if cds exon are successive (OK) or not (not valid)

    Parameters
    ---------- 
    cdsId:string
        cds Id
    blockList:list
    list of exon 

    Returns
    -------
    cut the program if the CDS exons are not successive  

    '''
    for i in range (0,len(blockList)-1):
		if blockList[i][1] !=  blockList[i+1][0]:
			print cdsId, ' invalid splicing structure (exon do not cover whole CDS): ', blockList, '\n'
			exit(-1)

def extractIntronGene(geneId, cdsId, blockList):
    """
    This function extracts introns list of gene based on a cds list of exons.

    Parameters
    ----------

    geneId: string
        gene Id
    cdsId:string
        cds Id
    blockList: list
        list of blocks


    Returns
    -------
    intronList:list
        list of intron
    """


    exonList = ''
    coord = []
    intronList = []
    for block in blockList:
        exonList = block[2:]
        coord.append(exonList[0])
        coord.append(exonList[1])

    for i in range(1, len(coord) - 2, 2):
        beginIntron = coord[i]
        endIntron = coord[i + 1]
        intronList.append([beginIntron, endIntron])

    return intronList



def poolCompare(gene_cds, cdsExon,  sourceData, cds2GeneExon,intronList,\
                    geneExon,step, compareExon,pairwiseMethod):
    """
    This function launchs comparison according to the method choosed to do the alignment and returns analysed results

    Parameters
    ----------

    gene_cds: list
            list of all pairs of gene and CDS information (id, sequence, ...)
    cdsExon: dictionary
        dictionay with CDS id as key and list of exons  as value
    sourceData: list
    	list that contains information about CDS : CDS ID, CDS sequence, gene ID, and exon list
    cds2GeneExon: dictionary
        dictionay with cds id as key and exons blocks list as values
    intronList:dictionary
               dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene
    geneExon: list
        list of exon of gene

    step: 

    Returns
    -------
    q: list
        list of information about the alignment between cds and gene
    """
    
    gene = gene_cds[0]
    cds = gene_cds[1]
    cdsId, cdsSeq, cdsGeneId, blockList_init = cds
    cdsLen = len(cdsSeq)

    geneId, geneSeq = gene
    geneLen = len(geneSeq)

    if(pairwiseMethod == "sfa"):
        #Local alignment
        blockList  = localAlignment(cdsId,cdsLen, cdsExon[cdsId],geneId, geneLen, geneSeq)
    
        #Block extension
        if  int(step) >= 2:
            blockList = Extend(cdsId, geneId, blockList, cdsLen, cdsSeq, cdsExon, geneLen, geneSeq, geneExon, intronList, compareExon)

        #Global alignment
        if  int(step) >= 3:
            blockList = AlignRestOfExons(cdsId, geneId, blockList, cdsExon, geneLen, cdsLen, geneExon, intronList, cdsSeq, geneSeq)
            
        
    if(pairwiseMethod == "splign"):
        cdsSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
        geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId  + '.fasta'
        splignOutputName = launch_splign(cdsSeqFile, geneSeqFile, cdsId, geneId )
        blockList = parse_splign_output_default(splignOutputName)

    if(cdsGeneId == geneId):
        blockList = blockList_init
    blockList = order(blockList)

        
    # analyse blocks results
    status, splicingSites, cdsTarget, texisting = analyse_result(blockList , cdsId, cdsSeq, geneId, geneSeq, sourceData, cds2GeneExon, cdsExon)
    
    # list to return
    q = [status, blockList, splicingSites, cdsTarget, texisting, geneId, cdsId]
    return q


    


def analyse_result(blockList, cdsId, cdsSequence, geneId, geneSequence, sourceData, cds2GeneExon, cdsExon):
    """
    This function analyses results of blocks found

    Parameters
    ----------


    blockList: list
    list of blocks found by the method of alignment
    geneId: string
            gene Id
    geneSequence: string
            gene sequence
    cds: string
        CDS Id
    cdsSequence: string
            CDS sequence

    cdsExon: dictionary
       dictionay with CDS id as key and exons interval list as value
    sourceData: list
       list that contain information about CDS like id CDS, CDS sequence and id of its gene

    cds2GeneExon: dictionary
       dictionay with cds id as key and exons blocks list as values


    Returns
    -------
    status: int

    splicing_sites: list
    list ofsplicing sites
    targetcds: string

    texisting: int

    """

    status = -1
    texisting = -1

   
    # create the target of the cds
    targetcds = create_target_cds(blockList, geneSequence)

    # compute the splicing sites
    splicing_sites = compute_splicing_sites(blockList, geneSequence)

    # create the target of exons
    targetexon = create_target_exon(blockList)

    
    if (len(blockList) > 0 and is_continuous_onCDS(blockList)):
        ccdsexon = cdsExon[cdsId]
        existing = False
        i = 0
        while(i < len(sourceData) and not existing):
            tcds = sourceData[i]
            tcdsid,tcdsseq,tgeneid,null = tcds
            ccdsexon = cdsExon[cdsId]
            tcdsexon = cdsExon[tcdsid]
            geneexon = cds2GeneExon[tcdsid]
            acceptor = True
            donor = True
            frame = True
            if(tgeneid == geneId and len(ccdsexon) == len(tcdsexon) == len(targetexon) == len(geneexon)):
                if(targetexon[0][1] != geneexon[0][1]):
                    donor = False
                for j in range(1,len(targetexon)-1):
                    if(targetexon[j][0] != geneexon[j][0]):
                        acceptor = False
                    if(targetexon[j][1] != geneexon[j][1]):
                        donor = False
                if(targetexon[-1][0] != geneexon[-1][0]):
                    acceptor = False
                for j in range(len(ccdsexon)):
                    sizeccds = ccdsexon[j][1]- ccdsexon[j][0]
                    sizetcds = tcdsexon[j][1]- tcdsexon[j][0]
                    if(sizeccds - sizetcds) % 3 != 0:
                        frame = False
            else:
                frame = False

            if(acceptor == True and donor == True and frame == True):
                existing = True
                status = STATUS_EXISTING_PROTEIN
                texisting = i
            i+=1
        if(not existing):
            if len(targetexon) == len(cdsExon[cdsId]):
                if(is_protein(targetcds)):
                    status = STATUS_PREDICTED_PROTEIN
                else:
                    status = STATUS_COMPLETE_CDS
            else:
                status = STATUS_COMPLETE_CDS_DIFFERENT_STRUCTURE
    else:
        status = STATUS_PARTIAL_CDS


    return status, splicing_sites,targetcds,texisting


def is_continuous_onCDS(blocklist):
	"""
	This function 

	Parameters
	----------

	blocklist:


	Returns
	-------
	boolean:
	"""
	for i in range(1,len(blocklist)):
		if(blocklist[i-1][1] < blocklist[i][0]):
		    return False
	return True



def create_target_cds(blocklist,geneseq):
    """
    This function 
    Parameters
    ----------

    blocklist:
    geneseq:
    

    Returns
    -------
    targetcds:
    
    """
    targetcds = ""
    
    for i in range(len(blocklist)):
        block = blocklist[i]
        targetcds +=  geneseq[block[2]:block[3]]

    return targetcds


def create_target_exon(blocklist):
    """
    This function 
    Parameters
    ----------

    blocklist:
        

    Returns
    -------
    targetexon:
    
    """
    targetexon = []

    for i in range(len(blocklist)):
        block = blocklist[i]
        targetexon.append(block[2:])

    return targetexon

def is_protein(targetcds):
    """
    This function return True if the predicted transcript have these 3 protein feature : length multiple of 3, start and stop codon
    Parameters
    ----------

    targetcds:
        

    Returns
    -------
    boolean:
    
    """

    length = len(targetcds)
        
    if(length%3 == 0 and targetcds[:3] == 'ATG' and (targetcds[length-3:] == 'TAA' or targetcds[length-3:] == 'TAG' or targetcds[length-3:] == 'TGA')):
        return True
    else:
        return False


def compute_splicing_sites(blocklist,geneseq):
	"""
	This function 
	Parameters
	----------

	blocklist:
	geneseq:


	Returns
	-------
	splicing_sites:

	"""
	splicing_sites = []
	for i in range(len(blocklist)-1):
		block1 = blocklist[i]
		block2 = blocklist[i+1]
		splicing_sites.append([geneseq[block1[3]:block1[3]+2],
			      geneseq[block2[2]-2:block2[2]]])
	return splicing_sites

def is_valid_splicing(splice):
    """
    This function 
    Parameters
    ----------

    
    splice:
        

    Returns
    -------
    
    
    """
    donor, acceptor = splice
    return ((donor == CANONICAL_DONOR and acceptor == CANONICAL_ACCEPTOR) or
            (donor == NC_DONOR1 and acceptor == NC_ACCEPTOR1) or
            (donor == NC_DONOR2 and acceptor == NC_ACCEPTOR2))
    
def all_valid_splicing_sites(splicing_sites):
    """
    This function 
    Parameters
    ----------

    
    splicing_sites:
        

    Returns
    -------
    all-valid
    
    """
    all_valid = True
    for splice in splicing_sites:
        all_valid = all_valid and is_valid_splicing(splice)
    return all_valid
