#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``step3_global`` **module description**:

This module performs the global alignment step of the method.

.. moduleauthor:: Safa Jammali



2017-2018

"""
from utils import *

MAXDPMATRIXSIZE= 10**7
GAP= -500
MATCH=  1000
MISMATCH= -500
Infinity = 10**7

minimumIntron =50
maximumIntron=2000

REAL_SPLICE_SITES =24000
ONE_REAL_SPLICE_SITES = 12000
REAL_EXON_JUNCTION = 24000

SPLICE_SIGNAL_GTAG = -10000
SPLICE_SIGNAL_GCAG = -15000
SPLICE_SIGNAL_ATAC = -20000
SPLICE_SIGNAL_OTHER = -50000

ThresholdIdentityMatch = 30

#########################
### alignement PD #######
#########################

def AlignRestOfExons(cdsId, geneId, blockList, cdsExon, geneLen, cdsLen, geneExon, intronList, cdsSeq, geneSeq ):
    """
    This function does the global alignment of remaining exons
    after steps 1 (local) and 2 (extend)

    Parameters
    ----------

    cdsId: string
        cds Id
    geneId: string
        gene Id
    blockListLocal: list
        current block list
    cdsExon :                
    geneLen: int
		gene length
    cdsLen: int
		cds length
    geneExon:
    intronList:
    cdsSeq: string
	sequence of CDS 
    geneSeq: string
	sequence of gene
   
    Returns
    -------
    blockList:list
        list of blocks
    """

    blockListDP = []
    blockListPDFound = []

    blockListDP = extractBlockListToDP(cdsId, geneId, blockList, cdsExon, geneLen, cdsLen)

    if  len(blockListDP) !=0:        
        blockList = InterBlocksPD(blockList, blockListDP, geneExon, cdsExon, intronList, cdsId, geneId, geneLen,
                                    cdsLen, cdsSeq, geneSeq)
   
    return blockList

def extractBlockListToDP(cdsId, geneId,blockListLocal, dictcdsExon,geneLen,cdsLen):
    """
    This function extract blocks (rest of exons that weren't aligned by the local alignement ) and that will be forced to be aligned

    Parameters
    ----------
    geneId: string
         gene id
     cdsId: string
         cds id

    blockListLocal: list
        list of blocks found by the local alignment
    dictcdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length
    Returns
    -------
    blockListPD:list
        list of blocks extracted and that will be extended


    """
    blockListPD = []
    listExonUnaligned = []
    cdsExon =     dictcdsExon [cdsId]

    if len(blockListLocal) == 0:
        blockListPD.append([0, cdsLen, 0, geneLen])
    else:
        for exon in cdsExon:
            unaligned = True
            for bloc in range(0, len(blockListLocal)):
                beginBlocCds = blockListLocal[bloc][0]
                endBlocCds = blockListLocal[bloc][1]
                if exon[1] > beginBlocCds and endBlocCds > exon[0]:
                    unaligned = False
            if(unaligned):
                listExonUnaligned.append(exon)
        for exon in listExonUnaligned:
            precb = [0, 0, 0, 0]
            nextb = [geneLen, geneLen, geneLen, geneLen]
            for bloc in range(0, len(blockListLocal)):
                if exon[0] >= blockListLocal[bloc][1]:
                    precb = blockListLocal[bloc]
                if exon[1] <= blockListLocal[bloc][0]:
                    nextb = blockListLocal[bloc]
                    break
            if (len(blockListPD) > 0 and blockListPD[-1][1] == exon[0]):
                blockListPD[-1][1] = exon[1]
                blockListPD[-1][3] = nextb[2]
            else:
                blockListPD.append([exon[0], exon[1], precb[3], nextb[2]])

    return blockListPD


    
def InterBlocksPD(resultsBloc,  blockListDP, geneExon,cdsExon, intronList,  cdsId, geneId, lenGene,lenCDS, cdsSeq, geneSeq ):
	"""
	This function extracts the inter blocks and launchs exact comparison for each interblock using dynamic prog.
	At the end, it recovers blocks and interblocks in the same list

	Parameters
	----------

	resultsBloc: list
			list contains current block lists
	geneExon: list
        list of exon of gene:
	cdsExon: list
		list of exons of cds
	intronList:dictionary
		dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene
	cdsId: string
		CDS Id
	geneId: string
		gene Id
	lenGene: int
		length of gene sequence
	lenCDS: int
		length of CDS sequence
	cdsSeq: string
		   CDS sequence
	lenGene:  string
			gene sequence


	Returns
	-------
	allBlocks: list
		list of bloks found for each pairs of CDS against gene
	"""

	
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
	posStartIntronGT=[]
	posEndIntronAG =[]

        for i in range(2):
            if(geneSeq[i: i+2] == "GT"):
                posStartIntronGT.append(i)
        for i in range(2,len(geneSeq)-2):
            if(geneSeq[i- 2: i] == "AG"):
                posEndIntronAG.append(i)
            if(geneSeq[i: i+2] == "GT"):
                posStartIntronGT.append(i)
        for i in range(len(geneSeq)-2,len(geneSeq)):
            if(geneSeq[i- 2: i] == "AG"):
                posEndIntronAG.append(i)
	
	extendedBlocks= order(blockListDP)
	# for each blocks
	
	for bloci in extendedBlocks:

            #identify extremity of two successifs blokcs
            [QueryStart, QueryEnd, SubjectStart, SubjectEnd]= bloci
		
		
            if (QueryEnd -QueryStart)* (SubjectEnd - SubjectStart) <= MAXDPMATRIXSIZE:
                if((SubjectEnd - SubjectStart) >0):
                    exonInGene, intronInGene, exonBlocksWithinCds,posExonStartGene, posExonEndGene, posExonStartCds, posExonEndCds, posStartIntron, posEndIntron \
                        = cdsGenerestrictionExonIntron(geneId, cdsId, geneExon, cdsExon, intronList, QueryStart, QueryEnd, \
												   SubjectStart, SubjectEnd)

                    

                    posStartIntronGTBloc=[]
                    for i in posStartIntronGT:
                        if(SubjectStart <= i <=SubjectEnd):
                            posStartIntronGTBloc.append(i)
                    posEndIntronAGBloc =[]
                    for i in posEndIntronAG:
                        if(SubjectStart <= i <=SubjectEnd):
                            posEndIntronAGBloc.append(i)
                    
                    newBloc = DynamicProgrammationAlignment(geneId, cdsId, geneSeq, cdsSeq,
                                                            QueryStart, QueryEnd, \
                                                            SubjectStart, SubjectEnd,\
                                                            geneExon, cdsExon, posExonStartGene,
                                                            posExonEndGene, posExonStartCds,
                                                            posExonEndCds, posStartIntron, posEndIntron, posStartIntronGTBloc, posEndIntronAGBloc, exonInGene, intronInGene, exonBlocksWithinCds)
                    
                    
                    #add the new block found
                    for  bloc in newBloc:
                        newBlocs += [bloc]


	# take all blocks
	allBlocks = resultsBloc + newBlocs

	#order blocks
	allBlocks= order(allBlocks)
	

	allBlocks = concatenate_exons(allBlocks, cdsSeq, geneSeq)
	
	return allBlocks


def cdsGenerestrictionExonIntron(geneId, cdsId,  geneExon, cdsExon, intronList, cdsBegin, cdsEnd, geneBegin, geneEnd):
	"""
	This function extracts the blocks of exons and introns of gene and exons of CDS belonging to the inter block

	Parameters
	----------

	geneId: string
		gene Id
	cdsId: string
		cds Id

	geneExon: list
		list of exon of gene
	cdsExon: list
		list of exons of cds

	intronList:dictionary
			   dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene

	cdsBegin: int
		value of the begining of the interblock at cds level
	cdsEnd: int
		end value  of the interblock at cds level
	geneBegin: int
		value of the beginig of the interblock at gene level
	geneEnd:int
		end value  of the interblock at gene level
	Returns
	-------

	posNucleotideInExonStart: dictionary
		dictionary of start positions of nucleotide in exon of gene belonging the interblock
	posNucleotideInExonEnd: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock
	posNucleotideInExonStartCds: dictionary
		dictionary of start positions of nucleotide in exon of cds  belonging the interblock
	posNucleotideInExonEndCds: dictionary
		dictionary of end positions of nucleotide in exon of cds  belonging the interblock
	posStartIntron: dictionary
		dictionary of start positions of nucleotide in intron  belonging the interblock
	posEndIntron:dictionary
		dictionary of end positions of nucleotide in intron  belonging the interblock

	exonInGene: dictionary
		 dictionary of list of exons in each gene, taking account the limit of interblock
	intronInGene: dictionary
		dictionary of list of introns in each gene, taking account the limit of interblock
	exonBlocksWithinCds: dictionary
		dictionary of list of exons in each cds, taking account the limit of interblock



	"""

	exonInGene={}
	posNucleotideInExonStart={}
	posNucleotideInExonEnd={}

	intronInGene = {}
	posStartIntron = {}
	posEndIntron = {}

	exonBlocksWithinCds = {}
	posNucleotideInExonStartCds = {}
	posNucleotideInExonEndCds = {}

	# extract interval block of exons  of gene belonging to the inter block
	exonInGene, posNucleotideInExonStart, posNucleotideInExonEnd =extracIntervalBlock(geneExon, geneId, exonInGene, \
												posNucleotideInExonStart, posNucleotideInExonEnd, geneBegin, geneEnd)

	# extract interval block of introns  of gene  belonging to the inter block
	intronInGene, posStartIntron , posEndIntron =extracIntervalBlock(intronList, geneId, intronInGene, posStartIntron,\
																	 posEndIntron, geneBegin, geneEnd)

	# extract interval block of exons  of CDS belonging to the inter block
	exonBlocksWithinCds, posNucleotideInExonStartCds, posNucleotideInExonEndCds, =extracIntervalBlock(cdsExon, cdsId, \
					exonBlocksWithinCds, posNucleotideInExonStartCds, posNucleotideInExonEndCds, cdsBegin, cdsEnd)


	return  exonInGene, intronInGene, exonBlocksWithinCds, posNucleotideInExonStart, posNucleotideInExonEnd, \
			posNucleotideInExonStartCds, posNucleotideInExonEndCds, posStartIntron, posEndIntron


def extracIntervalBlock(sourceDictionary, sourceId, blockInInterval, posNucleotideStartInInterval,\
						posNucleotideEndInInterval, beginInterval, endInterval):
	"""
	This function computes the blocks of exons and introns of gene and exons of CDS belonging to the inter block

	Parameters
	----------

	sourceDictionary: dictionary
		dictionary that contains list of cds exons or gene exons or gene introns
		 the keys are source Id and the values are the blocks interval of exon/intron source
	sourceId: string
		it can be cds id or gene Id

	blockInInterval: dictionary
		dictionary that contains list of cds exons or gene exons or gene introns beloning to the inter block
	posNucleotideStartInInterval: dictionary
		dictionary that cotains list of begin postion of cds exons or gene exons or gene introns beloning to the inter block
	posNucleotideEndInInterval: dictionary
		dictionary that cotains list of end postion of cds exons or gene exons or gene introns beloning to the inter block
	beginInterval: int
		begin of the inter block
	endInterval: int
		end of the inter block

	Returns
	-------
	blockInInterval:  dictionary
		 dictionary of list of exons/introns
	posNucleotideStartInInterval: dictionary
		dictionary of start positions of nucleotide
	posNucleotideEndInInterval: dictionary
		dictionary of end positions of nucleotide


	"""
	
        blockInInterval[sourceId] = []
        posNucleotideStartInInterval[sourceId] = []
        posNucleotideEndInInterval[sourceId] = []
	if(sourceId in sourceDictionary.keys()):
                values = sourceDictionary[sourceId]
                for i in range(0, len(values)):
                        minInter = values[i][0]
                        maxInter = values[i][1]
                        if (beginInterval <= minInter and minInter < endInterval):
                                if(beginInterval < maxInter and maxInter <= endInterval):
                                   blockInInterval[sourceId].append([minInter, maxInter])
                                   posNucleotideStartInInterval[sourceId].append(minInter)
                                   posNucleotideEndInInterval[sourceId].append(maxInter)

                                else:
                                   blockInInterval[sourceId].append([minInter, endInterval])
                                   posNucleotideStartInInterval[sourceId].append(minInter)
                        else:
                                if(beginInterval < maxInter and maxInter <= endInterval):
                                   blockInInterval[sourceId].append([beginInterval, maxInter])
                                   posNucleotideEndInInterval[sourceId].append(maxInter)
                                elif(minInter < beginInterval and maxInter > endInterval):
                                   blockInInterval[sourceId].append([beginInterval, endInterval])

	return blockInInterval, posNucleotideStartInInterval, posNucleotideEndInInterval


def concatenate_exons(blocklist,cdsseq,geneseq):
	"""
	This function concatenates all pairs of successive blocks that 
	whose cds_separation is less than MAX_CDS_SEPARATION_FOR_CONCAT, and 
	gene_separation doesnot exeed MAX_CDS_SEPARATION_FOR_CONCAT +/-MAX_SEPARATION_DIFF_FOR_CONCAT.
	Parameters
	----------

	blocklist:
	cdsseq:
	geneseq:

	Returns
	-------
	blocklist:
	"""


	block_to_delete=[]

	blocklist = order(blocklist)
    
	for i in range(len(blocklist)-1):

		block1_qs,block1_qe,block1_ss,block1_se=blocklist[i]
		block2_qs,block2_qe,block2_ss,block2_se=blocklist[i+1]

		cdsseparation = block2_qs-block1_qe
		geneseparation = block2_ss-block1_se

		if (cdsseparation <= MAX_CDS_SEPARATION_FOR_CONCAT):
		    
		    diff = abs(cdsseparation - geneseparation)
		    if (diff <= MAX_SEPARATION_DIFF_FOR_CONCAT and diff%3 == 0):
		        
		        alignment = pairwise2.align.globalms(geneseq[block1_se-5:block2_ss+5], cdsseq[block1_qe-5:block2_qs+5],2.0, 0.0, -5.0, -1.0)
		        alignementscore = -1 
		        for aln in alignment:
		            alignementscore = aln[4] 
		            
		        if (alignementscore > 0):
		            blocklist[i+1]=[block1_qs,block2_qe,block1_ss,block2_se]
		            block_to_delete.append(blocklist[i])

	for block in block_to_delete:
		blocklist.remove(block)

    	return blocklist





def DynamicProgrammationAlignment(geneId, cdsId, geneSequence, cdsSequence, \
								  cdsBeginSegment, cdsEndSegment, geneBeginSegment, geneEndSegment, \
								  geneExon, cdsExon, posExonStartGene, posExonEndGene, \
								  posExonStartCds, posExonEndCds, posStartIntron, posEndIntron, posStartIntronGT, posEndIntronAG, exonInGene, intronInGene,\
								  exonBlocksWithinCds):





	"""
	This function do a dynamic programmation alignment. It has  cases to compute and to choose the maximum one
	to fill it in the DP table

	Parameters
	----------

	cdsId: string
		CDS Id
	geneId: string
		gene Id
	cdsSequence: string
		   CDS sequence
	geneSequence:  string
		   gene sequence
	geneExon: dictionary
		dictionay with gene id as key and exons interval list as value
	cdsExon: dictionary
		dictionay with CDS id as key and exons interval list as value
	cdsBeginSegment: int
		value of begin segment in cds
	cdsEndSegment: int
		value of end segment in cds
	geneBeginSegment: int
		value of begin segment in gene
	geneEndSegment: int
		value of end segment in gene

	posExonStartGene: dictionary
		dictionary of start positions of nucleotide in exon of gene belonging the interblock
	posExonEndGene: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock
	posExonStartCds: dictionary
		dictionary of start positions of nucleotide in exon of cds  belonging the interblock
	posExonEndCds: dictionary
		dictionary of end positions of nucleotide in exon of cds  belonging the interblock
	posStartIntron: dictionary
		dictionary of start positions of nucleotide in intron  belonging the interblock
	posEndIntron:dictionary
		dictionary of end positions of nucleotide in intron  belonging the interblock

	exonInGene: dictionary
		 dictionary of list of exons in each gene, taking account the limit of interblock
	intronInGene: dictionary
		dictionary of list of introns in each gene, taking account the limit of interblock
	exonBlocksWithinCds: dictionary
		dictionary of list of exons in each cds, taking account the limit of interblock

	Returns
	-------
	newBlocs: list
		list of blocks exact aligned
	"""



	segmentCDS = '-' + cdsSequence[cdsBeginSegment: cdsEndSegment]
	segmentgene = '-' + geneSequence[geneBeginSegment: geneEndSegment]

	# create DP table
	M = create_matrix(len(segmentCDS), len(segmentgene))
	# create trace table
	Mtrace = create_matrix(len(segmentCDS), len(segmentgene))


	# initialization
	for i in range(1, len(segmentgene)):
		M[i][0] = i * GAP
		Mtrace[i][0] = [0, 0, 'casei']

	for j in range(1, len(segmentCDS)):
		M[0][j] = j * GAP
		Mtrace[0][j] = [0, 0, 'casej']

	# filling the DP table
	for j in range(1, len(segmentCDS)):
		for i in range(1, len(segmentgene)):

			i1 = 0
			j1 = 0
			iPrime = 0
			iMinusL = 0
			# calcul six cases in DP
			case1, case2, case3,  case5, iPrime, case6, iMinusL = case_alignment(geneSequence, \
																										  cdsSequence,
																										  M, i, j,
																										  geneExon,
																										  cdsExon,
																										  geneId, cdsId,
																										  posExonStartGene,
																										  posExonEndGene,
																										  posExonStartCds, \
																										  posExonEndCds,
																										  cdsBeginSegment,
																										  geneBeginSegment,
																										  segmentCDS,
																										  segmentgene, \
																										  posStartIntron,
																										  posEndIntron,
																										  posStartIntronGT,
																										  posEndIntronAG,
																										   exonInGene, \
																				intronInGene, exonBlocksWithinCds)




			# choose maximum case of the six
			M[i][j] = max(case1, case2, case3, case5, case6)
			# fill in the box according to the chosen case and store its trace
			if M[i][j] == case1:
				Mtrace[i][j] = [i - 1, j - 1, 'case1']
			elif M[i][j] == case2:
				Mtrace[i][j] = [i , j-1, 'case2']
			elif M[i][j] == case3:
				Mtrace[i][j] = [i-1, j , 'case3']
			elif M[i][j] == case5:
				Mtrace[i][j] = [iPrime, j, 'case5']
			elif M[i][j] == case6:
				Mtrace[i][j] = [iMinusL, j, 'case6']

	posmaxligne = foundPosition(M,segmentCDS, segmentgene )

	# display, identification and knowledge of the blocks alignment constructed by PD
	newBloc = traceAlignmentDP(posmaxligne, len(segmentCDS) - 1, Mtrace, geneId, cdsId, cdsBeginSegment,
							   geneBeginSegment, geneSequence, cdsSequence)

	# returns aligned blocks
	return newBloc


def create_matrix(m, n):
	"""


	Parameters
	----------

	m: int
		lengtn of cds
	n: int
		length of gene
	Returns
	-------
	M: list
		matrix n*m

	"""
	M = [[0 for i in range(m)] for j in range(n)]
	return M


def case_alignment(geneSequence, cdsSequence, M, i, j, geneExon, cdsExon, geneId, cdsId, posExonStartGene, posExonEndGene, \
				   posStartExonOnCds, posEndExonOnCds, cdsBeginSegment, geneBeginSegment, sequenceCdsInterblocks,\
				   sequenceGeneInterblocks, posStartIntron, posEndIntron, posStartIntronGT, posEndIntronAG, exonInGene, intronInGene, exonBlocksWithinCds):

	"""
	This function computes cases given in DP algorithms.

	Parameters
	----------

	cdsId: string
		CDS Id
	geneId: string
		gene Id
	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence

	geneExon: dictionary
		dictionay with gene id as key and exons interval list as value
	cdsExon: dictionary
		dictionay with CDS id as key and exons interval list as value

	cdsBeginSegment: int
		value of begin segment in cds
	geneBeginSegment: int
		value of begin segment in gene
	sequenceCdsInterblocks: string
		    CDS sequence in the interblock
	sequenceGeneInterblocks: string
		   origin CDS sequence in the interblock

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	posExonStartGene: dictionary
		dictionary of start positions of nucleotide in exon of gene belonging the interblock
	posExonEndGene: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock
	posExonStartCds: dictionary
		dictionary of start positions of nucleotide in exon of cds  belonging the interblock
	posExonEndCds: dictionary
		dictionary of end positions of nucleotide in exon of cds  belonging the interblock
	posStartIntron: dictionary
		dictionary of start positions of nucleotide in intron  belonging the interblock
	posEndIntron:dictionary
		dictionary of end positions of nucleotide in intron  belonging the interblock

	exonInGene: dictionary
		 dictionary of list of exons in each gene, taking account the limit of interblock
	intronInGene: dictionary
		dictionary of list of introns in each gene, taking account the limit of interblock
	exonBlocksWithinCds: dictionary
		dictionary of list of exons in each cds, taking account the limit of interblock

	Returns
	-------
	case1: int
		score of case1
	case2:int
		score of case2
	case3:int
		score of case3
	maxcas5: int
		score of case5
	iPrime: int
		position in gene, beging of real intron
	maxcas6:int
		score of case6
	iMinusL: int
		position in gene, beging of detected intron
	"""

	
	checkRealExonJunction = checkRealExonJunctionOnCDS(j+cdsBeginSegment-1, posStartExonOnCds[cdsId]+posEndExonOnCds[cdsId])


	
	# Compute case 5 score
	if  len(intronInGene[geneId]) == 0 :

		maxcas5 = -Infinity
		iPrime = 0

	else:
		maxcas5, iPrime = maxIntronscore(M, i, j, checkRealExonJunction, intronInGene[geneId], geneSequence, cdsSequence,\
						sequenceGeneInterblocks, sequenceCdsInterblocks, cdsBeginSegment, geneBeginSegment,\
						posExonStartGene[geneId], posExonEndGene[geneId], posStartIntron[geneId], posEndIntron[geneId])


	
	# Compute case 6 score
	maxcas6, iMinusL = maxIntronscore_case6(M, i, j, checkRealExonJunction, sequenceGeneInterblocks, sequenceCdsInterblocks,\
											geneBeginSegment, posStartIntron[geneId], posEndIntron[geneId],posStartIntronGT, posEndIntronAG)


	
	# Compute case 1 score
	case1 = case_1(M, i, j, MATCH, MISMATCH, sequenceGeneInterblocks, sequenceCdsInterblocks)
	# Compute case 2 score
	case2 = case_2(M, i, j, GAP)
	# Compute case 3 score
	case3 = case_3(M, i, j, GAP)
	
	#return scores
	return case1, case2, case3,  maxcas5, iPrime, maxcas6, iMinusL


def checkRealExonJunctionOnCDS(j,posStartExonOnCds ):

	"""
	This boolean function verif if actual postion in cds is a real junction or not.

	Parameters
	----------

	j: int
		actual position in CDS

	posExonStartCds: list
		list of start positions of nucleotide in exon of cds  belonging the interblock

	Returns
	-------
	boolean : True if yes ano False if not

	"""
	if j in posStartExonOnCds:
		 return True
	return False

def maxIntronscore(M, i, j, checkRealExonJunction, listIntronGene, geneSequence, cdsSequence, sequenceGeneInterblocks, \
				   sequenceCdsInterblocks, cdsBeginSegment, geneBeginSegment, posExonStartGene, posExonEndGene,\
				   posStartIntron, posEndIntron):


	"""
	This function computes score of case 5

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	checkRealExonJunction: boolean
		True if actual postion in cds is a real junction

	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence

	cdsBeginSegment: int
		value of begin segment in cds
	geneBeginSegment: int
		value of begin segment in gene
	sequenceCdsInterblocks: string
		    CDS sequence in the interblock
	sequenceGeneInterblocks: string
		   origin CDS sequence in the interblock

	posExonStartGene: dictionary
		dictionary of start positions of nucleotide in exon of gene  belonging the interblock

	posExonEndGene: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock

	posStartIntron: list
		list of start positions of nucleotide in exon of gene  belonging the interblock
	posEndIntron: list
		list of end positions of nucleotide in exon of gene  belonging the interblock


	listIntronGene: list
		 list of introns in each gene, taking account the limit of interblock

	Returns
	-------
	maxCsoreCase5: int
		score of case5
	iPrime: int
		position in gene, beging of real intron

	"""

	maxScoreCase5= -Infinity
	iPrime = 0

	scoreCase5 = -Infinity

	# compute score of real exon junction
	if checkRealExonJunction:
		scoreRealExonJunction = REAL_EXON_JUNCTION

	else:
		scoreRealExonJunction = 0

	# compute score of case 5 and its gene intron begining
	scoreSpliceSignals = SPLICE_SIGNAL_OTHER
	for intron in range(0, len(listIntronGene)):
		startIntron = listIntronGene[intron][0]
		endIntron = listIntronGene[intron][1]
		if i + geneBeginSegment -1 == endIntron  :

			scoreSpliceSites = checkRealSpliceSites(startIntron, endIntron, posStartIntron, posEndIntron)

			donor = sequenceGeneInterblocks[startIntron - geneBeginSegment +1: startIntron - geneBeginSegment + 2 +1]
			acceptor = sequenceGeneInterblocks[endIntron - geneBeginSegment - 2 +1: endIntron - geneBeginSegment +1]


			if donor== 'GT' and  acceptor == 'AG':
				scoreSpliceSignals = SPLICE_SIGNAL_GTAG
			elif donor == 'GC' and acceptor == 'AG':
				scoreSpliceSignals = SPLICE_SIGNAL_GCAG
			elif donor == 'AT' and acceptor == 'AC':
				scoreSpliceSignals = SPLICE_SIGNAL_ATAC

			scoreCase5 = M[startIntron - geneBeginSegment ][j] + scoreRealExonJunction + scoreSpliceSignals + scoreSpliceSites

			if maxScoreCase5 < scoreCase5:
				maxScoreCase5 = scoreCase5
				iPrime = startIntron - geneBeginSegment


		else:
			pass
		
	#return
	return maxScoreCase5, iPrime



def checkRealSpliceSites(i, l, posStartIntron, posEndIntron):

	"""
	This function verify if positions in gene are real start and real end intron in the gene
	and based on thes , it returns a score

	Parameters
	----------


	i: int
		it can be a start intron
	l: int
		it can be end of intron

	posStartIntron: list
		list of start positions of nucleotide in exon of gene  belonging the interblock
	posEndIntron: list
		list of end positions of nucleotide in exon of gene  belonging the interblock


	Returns
	-------
	scoreRealSpliceSites : int
		score of splice sites, it can be 0, 1000 or 2000

	"""
	scoreRealSpliceSites = 0
	if i in posStartIntron:
		scoreRealSpliceSites += ONE_REAL_SPLICE_SITES

	if l in  posEndIntron:
		scoreRealSpliceSites += ONE_REAL_SPLICE_SITES

	#return scores
	return scoreRealSpliceSites

#maxIntronscore_case6(M, i, j, checkRealExonJunction, sequenceGeneInterblocks, sequenceCdsInterblocks,\
#											geneBeginSegment, posStartIntron[geneId], posEndIntron[geneId])
#def maxIntronscore_case6(M, i, j, checkRealExonJunction, genesequence, CDSsequence, geneBeginSegment, posStartIntron, posEndIntron):
def maxIntronscore_case6(M, i, j, checkRealExonJunction, sequenceGeneInterblocks, sequenceCdsInterblocks, geneBeginSegment, posStartIntron, posEndIntron, posStartIntronGT, posEndIntronAG):
	"""
	This function computes score of case 6

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	checkRealExonJunction: boolean
		True if actual postion in cds is a real junction

	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence
	geneBeginSegment: int
		value of begin segment in gene



	posStartIntron: list
		list of start positions of nucleotide in exon of gene  belonging the interblock
	posEndIntron: list
		list of end positions of nucleotide in exon of gene  belonging the interblock


	Returns
	-------
	maxcas6: int
		score of case 6
	iMinusL: int
		position in gene, begining of intron


	"""

	maxcas6 = -Infinity
	iMinusL = 0
	scoreCase6 = -Infinity

	# compute score of real exon junction
	if checkRealExonJunction:
		scoreRealExonJunction = REAL_EXON_JUNCTION

	else:
		scoreRealExonJunction = 0

	# compute score of case 6 and its gene intron begining
	scoreSpliceSignals = SPLICE_SIGNAL_OTHER

        
        endIntron = i + geneBeginSegment -1
        acceptor = sequenceGeneInterblocks[i- 2: i]
        if(acceptor == 'AG'):
            # for l in range(minimumIntron,min(maximumIntron+1,i)):
	    #     startIntron = i + geneBeginSegment -1 -l
            #     donor = sequenceGeneInterblocks[i-l: i-l+2]
	    #     if donor== 'GT':
            for startIntron in posStartIntronGT:
                l = i + geneBeginSegment -1 - startIntron                
                if(l >= minimumIntron and l < min(maximumIntron+1,i)):
                    scoreSpliceSites = checkRealSpliceSites(startIntron, endIntron, posStartIntron, posEndIntron)                    
                    scoreSpliceSignals = SPLICE_SIGNAL_GTAG
                    scoreCase6 = M[i-1-l][j] + scoreRealExonJunction + scoreSpliceSignals + scoreSpliceSites
                    if maxcas6 < scoreCase6:
                        maxcas6 = scoreCase6
                        iMinusL = startIntron - geneBeginSegment
                        scoreSpliceSites = -Infinity

	return maxcas6, iMinusL



def case_1(M, i, j,  MATCH, MISMATCH, genesequence, CDSsequence):
	"""
	This function computes score of case 1

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	MATCH: int
		score of matching residu
	MISMATCH: int
		score of matching residu

	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence

	Returns
	-------
	res: int
		score of case 2
	"""


	if genesequence[i] == CDSsequence[j]:
		score = MATCH
	else:
		score = MISMATCH
	res = M[i - 1][j - 1] + score

	return res


def case_2(M, i, j, GAP):
	"""
	This function computes score of case 2

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	GAP: int
		score of gap


	Returns
	-------
	res: int
		score of case 2



	"""
	res = M[i ][j- 1] + GAP

	return res



def case_3(M, i, j, GAP):
	"""
	This function computes score of case 2

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	GAP: int
		score of gap


	Returns
	-------
	res: int
		score of case 3



	"""
	res = M[i- 1][j ] + GAP

	return res



def foundPosition(M,segmentCDS, segmentgene):
	"""
	This function computes the position that contains the max in the DP matrix M

	Parameters
	----------


	M: list
		DP matrix

	segmentCDS: string
		CDS sequence
	segmentgene: string
		gene sequence



	Returns
	-------
	posmaxligne: int
		position that contains the max
	"""

	# found max coord
	fincds = []
	



	for i in range(0, len(segmentgene)):

			fincds.append(M[i][len(segmentCDS) - 1])


	'''
	for j in range(0, len(segmentCDS)):

		for i in range(0, len(segmentgene)):

			if j == len(segmentCDS) - 1:
				fincds.append(M[i][j])
			#if i == len(segmentgene) - 1:
				#fingene.append(M[i][j])
	'''

	posmaxligne = 0

	maxfinc = -Infinity

	for j in range(0, len(fincds)):
		if fincds[j] > maxfinc:
			maxfinc = fincds[j]
                        posmaxligne = j
	
	return posmaxligne

def traceAlignmentDP(positionLine, positionColumn, M, geneId, cdsId, cdsBeginSegment, geneBeginSegment,geneSequence, cdsSequence):
	"""
	This function computes the list of blocks by doing the trace back on the DP matrix

	Parameters
	----------

	cdsId: string
		CDS Id
	geneId: string
		gene Id
	M: list
		DP matrix

	positionLine:

	positionColumn:


	cdsBeginSegment: int
		value of begin segment in cds

	geneBeginSegment: int
		value of begin segment in gene


	Returns
	-------
	newBlocs: list
		list of blocks exact aligned
	"""

        newBlocAccepted = []
	newBloc = [[0,0,0,0]]
	ancien =[]
	
	while (positionLine >= 0 or positionColumn >= 0):

		next =  M[positionLine][positionColumn]

                ancien = M[positionLine][positionColumn]

		[w, q, case] = M[positionLine][positionColumn]

		if w == 0 and q == 0:
                    if(case != 'case5' and case != 'case6'):
			newBloc[-1][0] = 0
                        newBloc[-1][2] = 0
                    else:
			newBloc[-1][0] = positionColumn
                        newBloc[-1][2] = positionLine
                    break
                if(case == 'case1'):
                    positionLine -= 1
                    positionColumn -=1
                if(case == 'case2'):
                    positionColumn -=1
                if(case == 'case3'):
                    positionLine -= 1
                if(case == 'case5' or case == 'case6'):
                    newBloc[-1][0] = positionColumn
                    newBloc[-1][2] = positionLine
                    newBloc.append([0,q,0,w+1])
                    positionLine = w
                    positionColumn = q

	for blocToTest in newBloc:
		[beginInCDS, endInCDS, beginInGene, endInGene] = blocToTest
                if(endInCDS != 0 and endInGene != 0):
                    if beginInCDS != 0:
			beginInCDS = beginInCDS - 1
                    if beginInGene != 0:
			beginInGene = beginInGene - 1
                    if endInCDS != 0:
			endInCDS = endInCDS - 1
                    if endInGene != 0:
			endInGene = endInGene - 1

                    sequence1 = geneSequence[beginInGene: endInGene]
                    sequence2 = cdsSequence[beginInCDS: endInCDS]
		
                    if len(sequence2)<=3:
			pass
                    else:

			alignment = launch_nwpairwise(sequence1, sequence2)
			
			identityMatch = computeAlignmentPercentIdentity(alignment[0][0], alignment[0][1])
			
			if identityMatch >= ThresholdIdentityMatch:
                            blocToTest = [beginInCDS+cdsBeginSegment,endInCDS+cdsBeginSegment,beginInGene+geneBeginSegment,endInGene+geneBeginSegment]
                            newBlocAccepted += [blocToTest]
	return newBlocAccepted

        


def launch_nwpairwise(sequence1, sequence2):
    """
    This function launches global pairwise alignment
    Parameters
    ----------

    sequence1:
    sequence2:
    

    Returns
    -------
    alignment:
   
    """
    alignment = pairwise2.align.globalms(sequence1, sequence2,2,0,-5,-1)
    return alignment

