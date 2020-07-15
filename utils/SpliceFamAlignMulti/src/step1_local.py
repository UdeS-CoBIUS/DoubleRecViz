#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``step1_local.py`` **module description**:

This module performs the local alignment step of the method

.. moduleauthor:: Safa Jammali



2017-2018

"""
from utils import *
from numpy import *

MAX_EVALUE = 1.0/10000000 # 10-7
LIMIT_ADD_NUCLEOTIDES = 3
INTERVAL_BLOC_COMPATIBLE_GENE_CDS = 40


###############################
### LOCAL COMPARISON ##########
###############################

def localAlignment(cdsId,cdsLen, cdsExon,geneId, geneLen, geneSeq):
	"""
	This function does the local alignment between CDS and gene 
        and returns blocks, that are significant hits between CDS and gene

	Parameters
	----------

	cdsId: string
		cds id
	cdsLen: int
		cds length
	cdsExon: list
		list of cds exons
	geneId: string
		gene id
	geneLen: int
		gene length
	geneSeq: string
		gene sequence

    Returns
      -------
      blockList: list
          lists of alignment blocks 
      """

	blockList = []

	evalue = MAX_EVALUE

	# recover cds and gene files name
	cdsFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
	geneFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId + '.fasta'

	# calls tblastx to do alignment
	blockList = launchAndTrimTblastx(cdsFile, geneFile, evalue, cdsId, geneId,cdsLen, geneLen, cdsExon, geneSeq)

        while(len(blockList) == 0 and evalue < 0.1):
                evalue = evalue*10
                blockList = launchAndTrimTblastx(cdsFile, geneFile, evalue, cdsId, geneId,cdsLen, geneLen, cdsExon, geneSeq)
	if (len(blockList) == 0):
		return []

	# returns list of hits if found, and returns [] if not found
	return blockList

def launchAndTrimTblastx( cdsFile, geneFile, evalue, cdsId, geneId,cdsLen, geneLen, cdsExon, geneSeq):
	"""
	This function calls TBlastX recursively to do the local alignment 
        between CDS and gene and returns blocks

	Parameters
	----------

	cdsFile: file name
		cds fasta file name 
	geneFile: file name
		gene fata file name
	evalue: float
		cutt off value for tblastx significant result
	cdsId: string
		cds id
	geneId: string
		gene id
	cdsLen:int
		length of cds
	geneLen:int
		length of gene
	cdsExon: list
		list of cds exons
	geneSeq: string
		gene sequence

    Returns
	-------
	blockList: list
		lists of alignment blocks 
	"""


        
	blastOutput = launch_tblastx(cdsFile, geneFile, evalue, [0, 0, 0, 0],[cdsLen, cdsLen, geneLen, geneLen], cdsId, geneId)

	blockList = parseTblastxTrimBlock(cdsId, geneId, cdsExon, blastOutput, geneSeq,evalue)
	os.remove(blastOutput)
	blockList = order(blockList)

	# returns list of hits
	return blockList

def parseTblastxTrimBlock(cdsId, geneId, cdsExon, blastOutput, geneSeq,evalue):

	"""
	This function takes tblastx results, check hits compatibility 
        and returns blocks that are compatible

	Parameters
	----------

	cdsId: string
		cds id
	geneId: string
		gene id
	cdsExon: list
		list of cds exons
	blastOutput: filename
		filename containing tblastx results (hits informations)
	blockQueryEnd : int
			value of the end residu to adjust blast result on the query
	blockSubjectEnd: int
			value of the end residu to adjust blast result on the subject
	geneSeq: string
		gene sequence
	evalue: float
		cutt off value for tblastx significant result

	Returns
	-------
	blocks: list
		lists of alignment blocks 
	"""

	
	dictExonBlocks = {}
	dictExonBlocksToaffect = {}
	incompatibledictExonBlocksToaffect = {}
	compatibledictExonBlocksToaffect = {}
	blocks = []
	valueExonBlocksNew = []

	
	firstExon= cdsExon[0]
	lastExon = cdsExon[len(cdsExon)-1]

	listExonOfCDS =[]
	#create an empty dictionary for each exon of cds
	for bloc in cdsExon:
		listExonOfCDS.append([bloc[0],bloc[1]])
		dictExonBlocks[(bloc[0],bloc[1])] =  []
		dictExonBlocksToaffect[(bloc[0],bloc[1])] =  []
		incompatibledictExonBlocksToaffect[(bloc[0],bloc[1])] =  []
		compatibledictExonBlocksToaffect[(bloc[0], bloc[1])] = []
	listHits=[]
        listHits= readBlastOut(blastOutput, 0, 0,evalue)

	listHits = order(listHits)
	
	for hit in listHits:
				# construct a hit

				beginNewHitCDS= hit[0]
				endNewHitCDS = hit[1]
				beginNewHitgene = hit[2]
				endNewHitgene = hit[3]
				evalueOfNewHit = hit[4]
	
				#affect hit in dictExonBlocks

				for indiceExonOfCDS in range (0, len(listExonOfCDS)):
					afterexonstart = beginNewHitCDS >= listExonOfCDS[indiceExonOfCDS][0]
					beforeexonend = beginNewHitCDS < listExonOfCDS[indiceExonOfCDS][1]

                                        #limit to remaining exons since hits are ordered
					if afterexonstart and beforeexonend:
						listExonOfCDS = listExonOfCDS[indiceExonOfCDS:len(listExonOfCDS)]
						break

	
				#  determine exon that contains the new hit
				maxLengthIntersection = -1
				segmentIntersection= None
				beginExonChoice = -1
				endExonChoice = -1
				i = 0
				continue_ = True
				
				while (continue_ == True and i <= len(listExonOfCDS)-1):
					e = listExonOfCDS[i]

					[beginExon, endExon] = e
					if beginExon >= endNewHitCDS:
						continue_ = False
						break



					segmentIntersectionLast, lengthIntersection = determineIntersection(beginExon, endExon,\
						int(beginNewHitCDS), int(endNewHitCDS), int(beginNewHitgene), int(endNewHitgene),cdsId, geneId,\
																					   geneSeq, firstExon, lastExon)


					if lengthIntersection > maxLengthIntersection:
							maxLengthIntersection = lengthIntersection
							segmentIntersection = segmentIntersectionLast

							[beginExonChoice , endExonChoice]= e
					i = i+1
				
				dictExonBlocks[(beginExonChoice, endExonChoice)].append((segmentIntersection, evalueOfNewHit))
	
	for exontoTreat in cdsExon:
		
		[beginExon, endExon] = exontoTreat
		listHitsExon = dictExonBlocks[(beginExon, endExon)]
		
		if len(listHitsExon) == 1 :			
			[(NewHit, evalueOfNewHit)] = listHitsExon
			dictExonBlocksToaffect[(beginExon, endExon)].append((NewHit, evalueOfNewHit))
		elif(len(listHitsExon) > 1):
			listHitsExonInOrder =sorted(listHitsExon, key=lambda colonnes: colonnes[1])
			for hitExon in listHitsExonInOrder:
                                (NewHit, evalueOfNewHit) = hitExon
				if len(dictExonBlocksToaffect[(beginExon, endExon)]) == 0:
					dictExonBlocksToaffect[(beginExon, endExon)].append((NewHit, evalueOfNewHit))
				else:
					dictExonBlocksToaffect[(beginExon, endExon)] = affectHit(dictExonBlocksToaffect[(beginExon, endExon)], NewHit, cdsId, geneId, evalueOfNewHit)

	
	for i in range(0, len(dictExonBlocksToaffect.keys())):
		
		for j in range(i + 1, len(dictExonBlocksToaffect.keys())):
			
			for hitexon in dictExonBlocksToaffect.values()[i]:
				exon1 = dictExonBlocksToaffect.keys()[i]
				for hitotherexon in dictExonBlocksToaffect.values()[j]:
					exon2 = dictExonBlocksToaffect.keys()[j]
					
					(hit1, eval1) = hitexon
					(hit2, eval2) = hitotherexon
					if exon1[1] <= exon2[0]:
						

						if int(hit1[3]) > int(hit2[2]):
							
							if (min(float(eval1), float(eval2)) == eval1):
								
								incompatibledictExonBlocksToaffect[exon2].append(hitotherexon)  #
							else:
								
								incompatibledictExonBlocksToaffect[exon1].append(hitexon)
					if exon2[1] <= exon1[0]:
						

						if int(hit2[3]) > int(hit1[2]):
							
							if (min(float(eval1), float(eval2)) == eval1):
								incompatibledictExonBlocksToaffect[exon2].append(
									hitotherexon)  
							else:
								incompatibledictExonBlocksToaffect[exon1].append(
									hitexon) 
	for cle, val in incompatibledictExonBlocksToaffect.items():
		
		if len(val) == 0:
			
			if len(dictExonBlocksToaffect[cle]) != 0:
				for k in dictExonBlocksToaffect[cle]:
					compatibledictExonBlocksToaffect[cle].append(k)
		
		else:
			newVal = []
			for j in dictExonBlocksToaffect[cle]:
				newVal = dictExonBlocksToaffect[cle]

			for k in val:
				if k in newVal:
					newVal.remove(k)
			compatibledictExonBlocksToaffect[cle] = newVal
		

	for blockItem in compatibledictExonBlocksToaffect.values():
			
			if (len(blockItem) == 0):
				pass

			if(len (blockItem) == 1):


				[(hitToAdd, eval) ]= blockItem
				blocks.append(hitToAdd)
			else:

				for element in blockItem:

					(hitToAdd, eval) = element
					blocks.append(hitToAdd)
	
	return blocks



def determineIntersection(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, firstExon, lastExon):
	"""
	This function determine intersection between cds exon and hit, it returns the intersection with its length

	Parameters
	----------

	beginExon: int
		value of begin exon
	endExon: int
		value of end exon
	beginNewHit: int
		value of begin of hit in cds
	endNewHit: int
		value of end of hit in cds
	beginNewHitgene: int
		value of begin of hit in gene
	endNewHitgene: int
		value of end of hit in gene

	geneId: string
		gene id
	cdsId: string
		cds id


	Returns
	-------

	newhit: list
		list that contains the new information of intersection between hit and exon
	length: int
		value of the length of the intersection between hit and exon

	"""

	length = 0
	newhit = None

	if beginNewHit > endExon :
		pass

	elif beginExon > endNewHit:
		pass

	#if newHit in exon
	elif (beginExon <= beginNewHit and beginNewHit <= endExon) and (beginExon <= endNewHit and endNewHit <= endExon):


		verifDonor = verifyLengthAndSitesDonor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, lastExon)
		if ((endExon - endNewHit ) <= 	LIMIT_ADD_NUCLEOTIDES) or verifDonor:
			endInCDS = endExon
			endInGene = endNewHitgene + (endExon - endNewHit)

		else:
			endInCDS = endNewHit
			endInGene = endNewHitgene

		verifAcceptor = verifyLengthAndSitesAcceptor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, firstExon)
		if (abs(beginNewHit - beginExon) <= LIMIT_ADD_NUCLEOTIDES) or verifAcceptor:
			beginInCDS = beginExon
			beginInGene =  beginNewHitgene - (beginNewHit - beginExon)

		else:
			beginInCDS = beginNewHit
			beginInGene = beginNewHitgene

		newhit = [beginInCDS, endInCDS , beginInGene, endInGene]
		length = endInCDS - beginInCDS

		
	# if exon in newHit
	elif beginNewHit <= beginExon and beginExon <= endNewHit and beginNewHit <= endExon and endExon <= endNewHit:

		beginNewHitGeneLevel = beginNewHitgene + abs(beginExon - beginNewHit)
		endNewHitGeneLevel = endNewHitgene - abs(endNewHit - endExon)
		newhit = [beginExon, endExon, beginNewHitGeneLevel,endNewHitGeneLevel]
		length = endExon - beginExon


		
	# if beginExon <= beginNewHit <= endExon  < endNewHit
	elif ((beginExon <= beginNewHit and beginNewHit <= endExon) and endNewHit > endExon):

		endNewHitGeneLevel = endNewHitgene - abs(endNewHit - endExon)
		verifAcceptor = verifyLengthAndSitesAcceptor(beginExon, endExon, beginNewHit, endExon, beginNewHitgene,
													 endNewHitGeneLevel, cdsId, geneId, geneSeq, firstExon)
		if (abs(beginNewHit - beginExon) <= LIMIT_ADD_NUCLEOTIDES) or verifAcceptor :
			beginInCDS = beginExon
			beginInGene = beginNewHitgene - abs(beginNewHit - beginExon)

		else:
			beginInCDS = beginNewHit
			beginInGene = beginNewHitgene

		newhit = [beginInCDS, endExon, beginInGene, endNewHitGeneLevel]
		length = endExon - beginInCDS
		
	# if beginNewHit < beginExon  <= endNewHit  <= endExon
	elif (beginNewHit < beginExon and (beginExon <= endNewHit and endNewHit <= endExon)):
		beginNewHitGeneLevel = beginNewHitgene + abs(beginExon - beginNewHit)

		verifDonor = verifyLengthAndSitesDonor(beginExon, endExon, beginExon, endNewHit, beginNewHitGeneLevel,
											   endNewHitgene, cdsId, geneId, geneSeq, lastExon)


		if ((endExon - endNewHit) <= LIMIT_ADD_NUCLEOTIDES) or verifDonor:

			endInCDS= endExon
			endNewHitGeneLevel = endNewHitgene + abs(endExon - endNewHit)

		else:
			endInCDS = endNewHit
			endNewHitGeneLevel = endNewHitgene

		newhit = [beginExon, endInCDS, beginNewHitGeneLevel, endNewHitGeneLevel]
		length = endInCDS - beginExon
		

	else:
		
		print 'Discarded hit in  determineIntersection'\
			, beginExon,endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene
		exit(-1)
	return newhit, length


def verifyLengthAndSitesAcceptor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, firstExon):
	"""
    This function verify the length of hit , if it is larger than exon length and the acceptor is AG , it returns True if not it returns false

    Parameters
    ----------

    beginExon: int
        value of begin exon
    endExon: int
        value of end exon
    beginNewHit: int
        value of begin of hit in cds
    endNewHit: int
        value of end of hit in cds
    beginNewHitgene: int
        value of begin of hit in gene
    endNewHitgene: int
        value of end of hit in gene

    geneId: string
        gene id
    cdsId: string
        cds id

	geneSeq: string
    	gene sequence
    Returns
    -------

    boolean

    """

	lenExon =  endExon - beginExon
	lenHit = endNewHit - beginNewHit
	hitCoverExon = False
	if lenHit >= (lenExon/2):
		hitCoverExon = True



	if beginExon > beginNewHit:

		endAcceptor = beginNewHitgene + abs(beginNewHit - beginExon)
	elif beginNewHit > beginExon:
		endAcceptor = beginNewHitgene - abs(beginNewHit - beginExon)

	else:
		endAcceptor = beginNewHitgene
	beginAcceptor = endAcceptor -2

	acceptor = geneSeq[beginAcceptor: endAcceptor] #acceptor== 'AG' end intron
	if  ((hitCoverExon and acceptor == 'AG') or (hitCoverExon and ([beginExon, endExon] == firstExon))):
		return True
	else:
		return False



def verifyLengthAndSitesDonor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, lastExon):
	"""
    This function verify th length of hit if it is larger than exon length and the donor is GT, it returns True if not it returns false

    Parameters
    ----------

    beginExon: int
        value of begin exon
    endExon: int
        value of end exon
    beginNewHit: int
        value of begin of hit in cds
    endNewHit: int
        value of end of hit in cds
    beginNewHitgene: int
        value of begin of hit in gene
    endNewHitgene: int
        value of end of hit in gene

    geneId: string
        gene id
    cdsId: string
        cds id

	geneSeq: string
    	gene sequence
    Returns
    -------

    boolean

    """

	lenExon =  endExon - beginExon
	lenHit = endNewHit - beginNewHit

	hitCoverExon = False
	if lenHit >= (lenExon/2):
		hitCoverExon = True



	if endExon > endNewHit:
		beginDonor = abs(endExon - endNewHit) +endNewHitgene

	elif endNewHit > endExon:
		beginDonor = endNewHitgene - abs(endExon - endNewHit)
	else:
		beginDonor = endNewHitgene

	endDonor = beginDonor + 2


	donor = geneSeq[beginDonor: endDonor]  #donor== 'GT' start intron

	if ( (hitCoverExon and donor == 'GT') or (hitCoverExon and ([beginExon, endExon] == lastExon)) ):
		return True
	else:
		return False


def affectHit(dictExonBlocksNew, segmentIntersection, cdsId, geneId, evalueOfNewHit):
	"""
	This function  compares and treats compatibility and returns  the compatible hits

	Parameters
	----------

	dictExonBlocksNew: list
		list of compatibl hits
	segmentIntersection:list
		new hit that will be treat

	geneId: string
		gene id
	cdsId: string
		cds id
	evalueOfNewHit: int
		value of e value given by tblastx for the hit

	Returns
	-------

	dictExonBlocksNew: list
		 contains the update list of hits


	"""

	dictExonBlocksNew_ = []
	countHit = 0
	HitAndEvalue = 0
	
	while HitAndEvalue < len(dictExonBlocksNew):
		
		(currentHit, evalueHit) = dictExonBlocksNew[HitAndEvalue]

		
		minBeginCds = min(segmentIntersection[0], currentHit[0])
		minBeginGene = min(segmentIntersection[2], currentHit[2])

		maxBeginCds = max(segmentIntersection[0], currentHit[0])
		maxBeginGene = max(segmentIntersection[2], currentHit[2])

		maxEndCds = max(segmentIntersection[1], currentHit[1])
		maxEndGene = max(segmentIntersection[3], currentHit[3])

		minDistanceCDS = abs(minBeginCds - maxBeginCds)
		minDistanceGene = abs(minBeginGene - maxBeginGene)
		compatible = 2
		case = ''
		# currentHit is include in the newHit, verify compatibility and add newHit or not
		if segmentIntersection[0] <= currentHit[0] and currentHit[0] <= segmentIntersection[1] \
				and segmentIntersection[0] <= currentHit[1] and currentHit[1] <= segmentIntersection[1] \
				and segmentIntersection[2] <= currentHit[2] and currentHit[2] <= segmentIntersection[3] \
				and segmentIntersection[2] <= currentHit[3] and currentHit[3] <= segmentIntersection[3]:
			
			countHit = countHit + 1
			# verify compatibility and add newHit or not
			if abs(minDistanceCDS - minDistanceGene) <= INTERVAL_BLOC_COMPATIBLE_GENE_CDS:

				# add the newHit: segmentIntersection
				if countHit == len(dictExonBlocksNew):
					
					dictExonBlocksNew_.append((segmentIntersection, evalueOfNewHit))
					
					break
				else:
				
					HitAndEvalue = HitAndEvalue + 1
				
			else:
				
				if (min(evalueHit, evalueOfNewHit) == evalueHit):
					dictExonBlocksNew_.append((currentHit, evalueHit))
					break

				else:
					
					if countHit == len(dictExonBlocksNew):

						
						dictExonBlocksNew_.append((segmentIntersection, evalueOfNewHit))
						break
					# remove currentHit
					else:

						HitAndEvalue = HitAndEvalue + 1
					
		# currentHit is different to the newHit, add newHit
		elif (segmentIntersection[0] > currentHit[1] and segmentIntersection[2] > currentHit[3]) \
				or (currentHit[0] > segmentIntersection[1] and currentHit[2] > segmentIntersection[3]):
			

			countHit = countHit + 1
			if countHit == len(dictExonBlocksNew):
				
				dictExonBlocksNew_.append((segmentIntersection, evalueOfNewHit))
				dictExonBlocksNew_.append((currentHit, evalueHit))
				break
			else:
				dictExonBlocksNew_.append((currentHit, evalueHit))
				HitAndEvalue = HitAndEvalue + 1
			

		# currentHit and newHit are overlaping, verify compatibility and add newHit or not
		elif (minBeginCds == segmentIntersection[0] and minBeginGene == segmentIntersection[2]) \
				or (minBeginCds == currentHit[0] and minBeginGene == currentHit[2]):
			
			countHit = countHit + 1
			compatible = 1
		else:
			countHit = countHit + 1
			
			compatible = 0
			case = 'non compatibe'

		if compatible == 1:

			if abs(minDistanceCDS - minDistanceGene) <= INTERVAL_BLOC_COMPATIBLE_GENE_CDS:

				
				if countHit == len(dictExonBlocksNew):
					
					dictExonBlocksNew_.append(
						([minBeginCds, maxEndCds, minBeginGene, maxEndGene], min(evalueHit, evalueOfNewHit)))
					break
				else:
					
					segmentIntersection = [minBeginCds, maxEndCds, minBeginGene, maxEndGene]
					evalueOfNewHit = min(evalueHit, evalueOfNewHit)
					HitAndEvalue = HitAndEvalue + 1

				
			else:
				
				if (min(evalueHit, evalueOfNewHit) == evalueHit):
					dictExonBlocksNew_.append((currentHit, evalueHit))
					break
				else:
					

					if countHit == len(dictExonBlocksNew):
						
						dictExonBlocksNew_.append((segmentIntersection, evalueOfNewHit))
						break
					else:

						HitAndEvalue = HitAndEvalue + 1

		elif (compatible == 0 and case == 'non compatibe'):

			if (min(evalueHit, evalueOfNewHit) == evalueHit):
				dictExonBlocksNew_.append((currentHit, evalueHit))
				break
			
			else:
				
				if countHit == len(dictExonBlocksNew):
					
					dictExonBlocksNew_.append((segmentIntersection, evalueOfNewHit))
					break
				else:

					HitAndEvalue = HitAndEvalue + 1

				
	return dictExonBlocksNew_


