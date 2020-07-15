#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``step2_extend.py`` **module description**:

This module performs the block extension step of the method.

moduleauthor:: Safa Jammali
Université de Sherbrooke Canada

2017-2018

"""

from step3_global import extractBlockListToDP
from step3_global import cdsGenerestrictionExonIntron
from step3_global import launch_nwpairwise
from step3_global import concatenate_exons
from utils import *

CUTOFF_IDENTITY_EXTEND = 30

def Extend(cdsId, geneId, blockListLocal, cdsLen, cdsSeq, cdsExon, geneLen, geneSeq, geneExon, intronList, compareExon):
    """
    This function extend the extremities of exons segments that were not aligned by the local alignement

    Parameters
    ----------

    geneId: string
  		gene Id
    cdsId: string
		cds Id
    blockListLocal: list
        list of blocks found by the local alignment
    cdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length

    cdsSeq: string
	sequence of CDS 
    geneSeq: string
	sequence of gene

    geneExon:
    intronList:
   
    Returns
    -------
    blockList: list
        blocks from local alignment with blocks found by the extension

    """

    blockListextend = []
    blockListextended = []

    
    blockListextend = extractBlockListToextend(cdsId, geneId, blockListLocal, cdsExon[cdsId], geneLen, cdsLen)


    extendedExtremities = extendExtremities(cdsId, geneId, blockListextend, cdsSeq, geneSeq)

    blockList = concatenateExtremities(cdsId, geneId, blockListLocal, extendedExtremities)


    if(compareExon == "Yes"): 
        blockListDP = extractBlockListToDP(cdsId, geneId, blockList, cdsExon, geneLen, cdsLen)


        if  len(blockListDP) !=0:        
            blockList = findSimilarExonsPD(blockList,  blockListDP, geneExon,cdsExon, intronList,  cdsId, geneId, geneLen, cdsLen, cdsSeq, geneSeq )

    return blockList

def extractBlockListToextend(cdsId, geneId,blockListLocal, cdsExon,geneLen,cdsLen):
    """
    This function extract blocks (rest of exons segments that weren't aligned by the local alignement ) and that will be extend

    Parameters
    ----------

    geneId: string
		gene Id
    cdsId: string
		cds Id
    blockListLocal: list
        list of blocks found by the local alignment
    cdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length
    Returns
    -------
    blockListToextend:list
        list of blocks extracted and that will be extended
   
    """

    listOfBlock = blockListLocal
    blockListToextend = []
    for exon in cdsExon:

        for bloc in range(0, len(listOfBlock)):
            # if block does not intersect exon
            if exon[0] >= listOfBlock[bloc][1] or exon[1] <= listOfBlock[bloc][0]:
                pass
            # if block intersects exon 
            else:
                # if ends differ
                if exon[1] != listOfBlock[bloc][1]:
                    # if there is another block after
                    if ((bloc + 1) < len(listOfBlock)):
                        # if next block starts after exon end
                        if listOfBlock[bloc + 1][0] >= exon[1]:
                            blockListToextend.append(
                                [listOfBlock[bloc][1], exon[1], listOfBlock[bloc][3],
                                 listOfBlock[bloc + 1][2], 'endExon'])
                            #blockListToextend.append(
                                #[listOfBlock[bloc][1], listOfBlock[bloc + 1][0], listOfBlock[bloc][3],
                                 #listOfBlock[bloc + 1][2], 'endExon'])
                        # if next block starts before exon end
                        else:

                            blockListToextend.append(
                                [listOfBlock[bloc][1], listOfBlock[bloc + 1][0], listOfBlock[bloc][3],
                                 listOfBlock[bloc + 1][2], 'middle'])
                    # if there is no block after
                    else:
                        blockListToextend.append([listOfBlock[bloc][1], exon[1], listOfBlock[bloc][3], geneLen, 'endExon'])
                # if starts differ
                if exon[0] != listOfBlock[bloc][0]:
                    # if there is another block before
                    if bloc != 0:
                        # if preceding block ends before exon start
                        if listOfBlock[bloc - 1][1] <= exon[0]:
                            blockListToextend.append(
                                [exon[0], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                 listOfBlock[bloc][2], 'beginExon'])
                            #blockListToextend.append(
                                #[listOfBlock[bloc - 1][1], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                 #listOfBlock[bloc][2], 'beginExon'])
                        # if preceding block end before exon
                        else:
                            blockListToextend.append(
                                [listOfBlock[bloc - 1][1], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                 listOfBlock[bloc][2], 'middle'])
                    # if there is no block before
                    else:
                        blockListToextend.append([exon[0], listOfBlock[bloc][0], 0, listOfBlock[bloc][2], 'beginExon'])

    newblockListToextend = []
    for i in blockListToextend:
        if i not in newblockListToextend:
            newblockListToextend.append(i)
   
    return newblockListToextend


def extendExtremities(cdsId, geneId, blockListextend, cdsSeq, geneSeq):
    max_indel = 12
    extendedExtremities = []
    for bloc in blockListextend:
        [cdsBegin, cdsEnd, geneBegin, geneEnd, position] = bloc
        extend  = []
        cdsLength = cdsEnd - cdsBegin
        geneLength = geneEnd - geneBegin
        if(position == "middle"):
            if(abs(cdsLength - geneLength) <= 3*max_indel and (cdsLength - geneLength)%3 == 0):
                extend = [cdsBegin, cdsEnd, geneBegin, geneEnd, position]
                extendedExtremities.append(extend)
        elif(position == "endExon"):
            precId = CUTOFF_IDENTITY_EXTEND
            for i in range(min(max_indel,(geneEnd-geneBegin-cdsLength)/3)):
                cdsseqtemp = cdsSeq[cdsBegin:cdsEnd]
                geneseqtemp = geneSeq[geneBegin+3*i:geneBegin+cdsLength+3*i]
                splicesite = geneSeq[geneBegin+cdsLength+3*i:geneBegin+cdsLength+3*i+2]
                identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)
                if(splicesite == "GT" and identityAlignment > precId):
                    extend = [cdsBegin, cdsEnd, geneBegin+3*i, geneBegin+cdsLength+3*i, position, identityAlignment]
                    precId = identityAlignment
            for i in range(min(max_indel,cdsLength/3)):
                cdsseqtemp = cdsSeq[cdsBegin+3*i:cdsEnd]
                geneseqtemp = geneSeq[geneBegin:geneBegin+(cdsEnd-(cdsBegin+3*i))]
                splicesite = geneSeq[geneBegin+(cdsEnd-(cdsBegin+3*i)):geneBegin+(cdsEnd-(cdsBegin+3*i))+2]
                identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)
                if(splicesite == "GT" and identityAlignment > precId):
                    extend = [cdsBegin+3*i, cdsEnd, geneBegin, geneBegin+(cdsEnd-(cdsBegin+3*i)), position, identityAlignment]
                    precId = identityAlignment
            if(len(extend) != 0):
                extendedExtremities.append(extend)
        elif(position == "beginExon"):
            precId = CUTOFF_IDENTITY_EXTEND
            for i in range(min(max_indel,(geneEnd-geneBegin-cdsLength)/3)):
                cdsseqtemp = cdsSeq[cdsBegin:cdsEnd]
                geneseqtemp = geneSeq[geneEnd-cdsLength-3*i:geneEnd-3*i]
                splicesite = geneSeq[geneEnd-cdsLength-3*i-2:geneEnd-cdsLength-3*i]
                identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)
                if(splicesite == "AG" and identityAlignment > precId):
                    extend = [cdsBegin, cdsEnd, geneEnd-cdsLength-3*i, geneEnd-3*i, position, identityAlignment]
                    precId = identityAlignment
            for i in range(min(max_indel,cdsLength/3)):
                cdsseqtemp = cdsSeq[cdsBegin:cdsEnd-3*i]
                geneseqtemp = geneSeq[geneEnd-(cdsEnd-3*i-cdsBegin):geneEnd]
                splicesite = geneSeq[geneEnd-(cdsEnd-3*i-cdsBegin)-2:geneEnd-(cdsEnd-3*i-cdsBegin)]
                identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)
                if(splicesite == "AG" and identityAlignment > precId):
                    extend = [cdsBegin, cdsEnd-3*i, geneEnd-(cdsEnd-3*i-cdsBegin), geneEnd, position, identityAlignment]
                    precId = identityAlignment
            if(len(extend) != 0):
                extendedExtremities.append(extend)
            
    return extendedExtremities

def concatenateExtremities(cdsId, geneId, blockListLocal, extendedExtremities):
    blockList = []
    sortedlist = sorted(blockListLocal + extendedExtremities)
    i = 0
    while(i < len(sortedlist)):
        bloc = sortedlist[i]
        if(len(bloc) == 4):
            blockList.append(bloc)
        elif(bloc[4] == "middle"):
            blockList[-1] = [blockList[-1][0], sortedlist[i+1][1],blockList[-1][2],sortedlist[i+1][3]]
            i += 1
        elif(bloc[4] == "beginExon"):
            blockList.append([bloc[0],sortedlist[i+1][1],bloc[2],sortedlist[i+1][3]])
            i += 1
        elif(bloc[4] == "endExon"):
            blockList[-1] = [blockList[-1][0], bloc[1],blockList[-1][2],bloc[3]]
        i += 1
    return blockList

def findSimilarExonsPD(resultsBloc,  blockListDP, geneExon,cdsExon, intronList,  cdsId, geneId, lenGene,lenCDS, cdsSeq, geneSeq ):
    	newBlocs=[]

	exonInGene ={}
	intronInGene ={}
	exonBlocksWithinCds ={}
	posExonStartGene=[]
	posExonEndGene=[]
	posExonStartCds=[]
	posExonEndCds=[]
	posStartIntron=[]
	posEndIntron =[]

	
	extendedBlocks= order(blockListDP)
	# for each blocks
	
	for bloci in extendedBlocks:

		#identify extremity of two successifs blokcs
		[QueryStart, QueryEnd, SubjectStart, SubjectEnd]= bloci
		
                if((SubjectEnd - SubjectStart) >0):
                    exonInGene, intronInGene, exonBlocksWithinCds,posExonStartGene, posExonEndGene, posExonStartCds, posExonEndCds, posStartIntron, posEndIntron \
                        = cdsGenerestrictionExonIntron(geneId, cdsId, geneExon, cdsExon, intronList, QueryStart, QueryEnd, \
												   SubjectStart, SubjectEnd)


                    newBloc = findSimilarExons(exonBlocksWithinCds[cdsId],exonInGene[geneId], cdsSeq, geneSeq)

                    #add the new block found
                    for  bloc in newBloc:
                        newBlocs += [bloc]



	# take all blocks
	allBlocks = resultsBloc + newBlocs

	#order blocks
	allBlocks= order(allBlocks)
	

	allBlocks = concatenate_exons(allBlocks, cdsSeq, geneSeq)
	
	return allBlocks

def findSimilarExons(exonBlocksWithinCdsList,exonInGeneList, cdsSeq, geneSeq):

    newBlocks = []
    if(len(exonBlocksWithinCdsList) > 0 and len(exonInGeneList) > 0):
        identityMatrix = []
        for i in range(len(exonBlocksWithinCdsList)):
            identityMatrix.append([])
            exonCds = exonBlocksWithinCdsList[i]
            for j in range(len(exonInGeneList)):
                exonGene = exonInGeneList[j]
                sequence1 = geneSeq[exonGene[0]:exonGene[1]]
                sequence2 = cdsSeq[exonCds[0]:exonCds[1]]
                identityMatch = 0
                if len(sequence1)<=0 or  len(sequence2)<=0:
                    pass
                else:
                    alignment = launch_nwpairwise(sequence1, sequence2)
                    identityMatch = computeAlignmentPercentIdentity(alignment[0][0], alignment[0][1])
                identityMatrix[i].append(identityMatch)

        newBlocks = findMostSimilarExons(exonBlocksWithinCdsList,exonInGeneList, identityMatrix)
                        
    return newBlocks

def findMostSimilarExons(exonBlocksWithinCdsList,exonInGeneList, identityMatrix):
    newBlocks = []
    maxIdentity = 0
    iMatch = -1
    jMatch = -1
    if(len(identityMatrix) > 0 and len(identityMatrix[0]) > 0):
        for i in range(len(identityMatrix)):
            for j in range(len(identityMatrix[0])):
                if(identityMatrix[i][j] > maxIdentity):
                    maxIdentity = identityMatrix[i][j];
                    iMatch = i
                    jMatch = j
                    
    if maxIdentity >= CUTOFF_IDENTITY_EXTEND:
        newBlocks.append([exonBlocksWithinCdsList[iMatch][0],exonBlocksWithinCdsList[iMatch][1], exonInGeneList[jMatch][0],exonInGeneList[jMatch][1]])

        newBlocks += findMostSimilarExons(exonBlocksWithinCdsList[:iMatch],exonInGeneList[:jMatch], [line[:jMatch] for line in identityMatrix[:iMatch]]) + findMostSimilarExons(exonBlocksWithinCdsList[iMatch+1:],exonInGeneList[jMatch+1:], [line[jMatch+1:] for line in identityMatrix[iMatch+1:]])

    return newBlocks
