#!/usr/bin/python

'''
Authors : Esaie Kuitche 
Example : python

'''

import sys
import argparse
import os
import ntpath
import numpy as np
import pandas as pd
from math import sqrt
import random
import operator
from scipy.spatial import distance
import collections
import string
import copy
import itertools
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import pylab
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint
from Bio.Seq import Seq
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate

def chooseCminCmax(n, mappingGeneToTranscript):
	Cmin = 0
	for geneId, seqs in mappingGeneToTranscript.items():
		if len(seqs)>Cmin:
			Cmin = len(seqs)

	Cmax =  int(sqrt(n)) + 1
	if Cmax > Cmin:
		if Cmin == 1:
			Cmin = 2
		return Cmin, Cmax
	else:
		if Cmin == 1:
			Cmin = 2		
		return Cmax, Cmin

def initializeClusterCenters(clusters, X, x_coords):
	index = range(len(X))
	random.shuffle(index)	

	i = 0
	for key in clusters.keys():	
		#print(key, index[i])					
		initval = x_coords[X[index[i]]]
		#initval = [x*random.random() for x in initval]
		clusters[key] = initval
		i += 1
	#exit()
	return clusters

def compute_distance(x_1, x_2):
	#dist = distance.cosine(x_1, x_2)
	dist = distance.euclidean(x_1, x_2)
	return dist

def FCM(U, fuzzifier, X, x_coords, clusterCentersIndex, clusterCenters, epsilon):			
	#print(U)
	while True:
		#print("______________________________________________________")
		for x in X:
			x_k = x_coords[x]
			for v in clusterCentersIndex:
				v_i = clusterCenters[v]
				#print(x_k, v_i, compute_distance(x_k, v_i))
				if compute_distance(x_k, v_i) == 0:
					u_k_i = 1
				else:
					flag = False
					for j in clusterCentersIndex:

						if j == x :
							pass
						else:
							v_j = clusterCenters[j]
							if compute_distance(x_k, v_j)==0:							
								u_k_i = 0
								flag = True
								break
					if flag == True:
						pass
					else:
						u_k_i = 0.0
						for j in clusterCentersIndex:
							v_j = clusterCenters[j]
							##print(compute_distance(x_k, v_i)/compute_distance(x_k, v_j), pow((compute_distance(x_k, v_i)/compute_distance(x_k, v_j)), 2.0/(fuzzifier-1)))
							#print(x_k, v_i, compute_distance(x_k, v_i), x_k, v_j, compute_distance(x_k, v_j))
							u_k_i += pow((compute_distance(x_k, v_i)/compute_distance(x_k, v_j)), 2.0/(fuzzifier-1))
						#print("--------------------------------------", u_k_i)
						u_k_i = 1.0/u_k_i					
				#print(u_k_i)
				U[X.index(x), clusterCentersIndex.index(v)] = u_k_i
		#print(U)
		clusterCenters_1 = {}

		for v in clusterCentersIndex:
			num = 0.0
			den = 0.0
			for x in X:
				x_k = x_coords[x]
				num += np.dot(pow(U[X.index(x), clusterCentersIndex.index(v)], fuzzifier),x_k)
				den += pow(U[X.index(x), clusterCentersIndex.index(v)], fuzzifier) + 0.0001
			clusterCenters_1[v] = num/den

		flag = True
		for v in clusterCentersIndex:
			ratio = compute_distance(clusterCenters[v], clusterCenters_1[v])/np.linalg.norm(clusterCenters_1[v])
			if ratio >= epsilon:
				flag = False
				break		
		if flag == False:			
			clusterCenters = clusterCenters_1			
		else:
			break
	#print(U)
	return clusterCenters, U

def extractCluster(U, clusterCenters, clusterCentersIndex, X, mappingGeneToTranscript):
	mappingGeneToTranscript = detailGeneTranscript(mappingGeneToTranscript)	
	clusters = {}
	for clusterCenter in clusterCentersIndex:
		clusters[clusterCenter] = []

	
	allMenberShip = {}
	for x in X:
		x_k = X.index(x)
		max = {}
		indexClusterMax = 0
		for v in clusterCentersIndex:
			v_i = clusterCentersIndex.index(v)
			max[v] = U[x_k,v_i]				
		max = sorted(max.items(), key=operator.itemgetter(1), reverse=True)

		e = max[0]
		cluterCenter = e[0]
		menbership = e[1]
		clusters[cluterCenter].append(x)
		allMenberShip[x] = max
	#print(clusters)
	for cluster, seq in clusters.items():	
		seq_final = []
		genes = [mappingGeneToTranscript[x] for x in seq]
		repeatedGenes = [item for item, count in collections.Counter(genes).items() if count > 1]		
		if len (repeatedGenes) > 0:
			clusterIdOfX = {}
			seqAll = copy.deepcopy(seq)
			for repeatedGene in repeatedGenes:
				ind = [i for i, j in enumerate(genes) if j in repeatedGene]	
				#print (ind, seq)			
				conflicts_x = [seqAll[x] for x in ind]
				maxMenbership = 0.0
				for x_i in conflicts_x:			
					index = -1
					allMenberShipTmp = allMenberShip[x_i]

					for e in allMenberShipTmp:
						if cluster == e[0]:						
							index = allMenberShipTmp.index(e)
							clusterIdOfX[x_i] = index
							break
					menbershipValues = allMenberShipTmp[index][1]

					if menbershipValues > maxMenbership:
						maxMenbership = menbershipValues
						x_max = x_i

					
				remaing_conflicts_x = [x for x in conflicts_x if x != x_max]

				seq = [x for x in seq if x not in remaing_conflicts_x]

				for x_i in remaing_conflicts_x:

					for i in range(clusterIdOfX[x_i] +1, len(clusterCenters)):
						nextClass = clusters[allMenberShip[x_i][i][0]]

						nextClassgenes = [mappingGeneToTranscript[x] for x in nextClass]

						currentGene = mappingGeneToTranscript[x_i]
						if currentGene in nextClassgenes:
							pass
						else:
							nextClass.append(x_i)
							break
		clusters[cluster] = seq
	
	#print(U)
	return clusters


def detailGeneTranscript(mappingGeneToTranscriptInit):
	mappingGeneToTranscript = {}
	for k, v in mappingGeneToTranscriptInit.items():
		for seqId in v:
			mappingGeneToTranscript[seqId] = k
	return mappingGeneToTranscript


def validityIndex(X, x_coords, clusterCentersIndex, clusterCenters, U):
	n = len(X)
	nb_cluster = len(clusterCentersIndex)
	all_x = []
	for x in X:
		all_x.append(x_coords[x])
	x_mean = np.mean(all_x, axis=0)
	x_sigma = []

	for p in range(len(x_mean)):
		x_sigma_tmp = 0.0
		for k in range(n):
			x_sigma_tmp += ((x_coords[X[k]][p] - x_mean[p])*(x_coords[X[k]][p] - x_mean[p]))
		x_sigma_tmp /= n
		x_sigma.append(x_sigma_tmp)

	v_sigma = {}

	for i in range(nb_cluster):		
		v_i = clusterCentersIndex[i]
		v_sigma[v_i] = []
		v_sigma_tmp_list = []
		for p in range(len(x_mean)):
			v_sigma_tmp = 0.0
			for k in range(n):
				u_k_i = U[k,i]
				v_sigma_tmp +=  u_k_i*((x_coords[X[k]][p] - clusterCenters[v_i][p])*(x_coords[X[k]][p] - clusterCenters[v_i][p]))				
			v_sigma_tmp /= n
			v_sigma_tmp_list.append(v_sigma_tmp)	
		v_sigma[v_i] = v_sigma_tmp_list

	scat_c = 0.0
	scat_c_num = 0.0
	scat_c_den = 0.0

	for v_i, vect in v_sigma.items():
		scat_c_num += np.linalg.norm(vect)

	scat_c_num /= nb_cluster

	scat_c_den = np.linalg.norm(x_sigma)

	scat_c = scat_c_num/scat_c_den

	sep_c = 0.0

	d_min = 999999999.0
	d_max = 0.0

	#print (clusterCenters)

	for i in range(nb_cluster):
		v_i = clusterCenters[clusterCentersIndex[i]]
		tmp_sum = 0.00000000001
		for j in range(nb_cluster):
			v_j = clusterCenters[clusterCentersIndex[j]]
			d = compute_distance(v_i, v_j)			
			#print (i, j, v_i, v_j, d)
			if i != j:			
				tmp_sum += d
				if d < d_min:
					d_min = d
				if d > d_max:
					d_max = d
		sep_c += 1/tmp_sum

	#print (d_min, d_max)
	d_min += 0.00000000001

	sep_c = ((d_max*d_max)/(d_min*d_min))*sep_c	

	validity_index = scat_c + sep_c

	return scat_c, sep_c

def findWorstClusterToSplit(U, clusters, X, x_coords, clusterCentersIndex, clusterCenters):

	min = 999999.0
	min_cluster = []

	for clusterCenter, sequences in clusters.items():
		#print(clusterCenter, sequences)
		if len(sequences)<2:
			pass
		else:
			clusterScore = 0.0
			k = clusterCentersIndex.index(clusterCenter)
			for seq in sequences:
				i = X.index(seq)
				clusterScore += U[i, k]

			clusterScore /= len(sequences)
			if clusterScore < min :
				min = clusterScore
				min_cluster = [clusterCenter, sequences]
	if min_cluster == []:
		return []
	else:
		return min_cluster[0], min_cluster[1]


def findNextTwoCenters(clusterCenterToSplit, sequencesToSplit, clusterCentersIndex, clusterCenters, U, x_coords):
	char_set = string.ascii_uppercase + string.digits

	v_io = clusterCenterToSplit
	E = sequencesToSplit
	bannedSeq = []
	maxDistance = []
	max	= 0.0
	v_i1 = ""
	v_i2 = ""
	#print(sequencesToSplit)
	for v_i1_tmp in sequencesToSplit:
		E0 = []
		E1 = []
		if v_i1_tmp in bannedSeq:
			pass
		else:			

			distanceToOtherClusterCenters = 0.0
			for clusterCenter in clusterCentersIndex:
				distanceToOtherClusterCenters += compute_distance(x_coords[v_i1_tmp], clusterCenters[clusterCenter])

			if distanceToOtherClusterCenters >= max:
				for seq in sequencesToSplit:					
					#print(x_coords[v_i1_tmp], x_coords[seq], clusterCenters[v_io], x_coords[seq], compute_distance(x_coords[v_i1_tmp], x_coords[seq]), compute_distance(clusterCenters[v_io], x_coords[seq]))
					if compute_distance(x_coords[v_i1_tmp], x_coords[seq]) <= compute_distance(clusterCenters[v_io], x_coords[seq]):
						E1.append(seq)
					else:
						E0.append(seq)
				#print(E1, E0)
				if (len(E1)*1.0)/len(E) >= 0.1:
					max = distanceToOtherClusterCenters
					v_i1 = v_i1_tmp
	bannedSeq.append(v_i1)

	max	= 0.0
	#print(sequencesToSplit, bannedSeq)

	for v_i2_tmp in sequencesToSplit:		

		E1 = []
		E2 = []
		if v_i2_tmp in bannedSeq:
			pass
		else:						
			distanceToOtherClusterCenters = 0.0
			for clusterCenter in clusterCentersIndex:
				distanceToOtherClusterCenters += compute_distance(x_coords[v_i2_tmp], clusterCenters[clusterCenter])
			
			if distanceToOtherClusterCenters >= max:	
							
				for seq in sequencesToSplit:					
					if compute_distance(x_coords[v_i2_tmp], x_coords[seq]) <= compute_distance(x_coords[v_i1], x_coords[seq]):
						E2.append(seq)
					else:
						E1.append(seq)

				if (len(E2)*1.0)/len(E) >= 0.1:
					max = distanceToOtherClusterCenters
					v_i2 = v_i2_tmp
					bannedSeq.append(v_i2_tmp)


	#print(v_i1, v_i2, x_coords)
	coord_v_i1 = x_coords[v_i1]
	coord_v_i2 = x_coords[v_i2]

	v_i1 = v_i1 + "_" + ''.join(random.sample(char_set*6, 6)) 
	v_i2 = v_i2 + "_" + ''.join(random.sample(char_set*6, 6))

	return v_i1, coord_v_i1, v_i2, coord_v_i2, E1, E2

	
def selectAndCompleteTranscriptTree(finalClusters, finalClusterCentres, finalclusterCentersIndex, mappingGeneToTranscript, X, x_coords, U):
	allGenes = mappingGeneToTranscript.keys()
	mappingGeneToTranscriptDetail = detailGeneTranscript(mappingGeneToTranscript)	

	maxElement = []
	maxSize = 0
	for clusterId, seqs in finalClusters.items():
		#print(clusterId, seqs)
		if len(seqs) > maxSize:
			maxElement = [[clusterId, seqs]]
			maxSize = len(seqs)
		elif len(seqs) == maxSize:
			maxElement.append([clusterId, seqs])
	
	finalClustersComplete = []
	for cluster in maxElement:
		center = cluster[0]
		seqs = cluster[1]
		genesOfSeq = [mappingGeneToTranscriptDetail[x] for x in seqs]
		missingGenes = [x for x in allGenes if x not in genesOfSeq]		
		seqIdToAdd = []		

		for geneId in missingGenes:
			
			potentialMissingSeqID = mappingGeneToTranscript[geneId]
			minDistance = 99999999.0
			maxMemberShipIds = []
			for seqId in potentialMissingSeqID:
				i = X.index(seqId)
				
				k = finalclusterCentersIndex.index(center)
				dist = compute_distance(x_coords[seqId], finalClusterCentres[center])
				#print (seqId, center, dist)				

				if dist < minDistance:
					maxMemberShipIds = [seqId]
					minDistance = dist
				elif dist == minDistance:
					maxMemberShipIds.append(seqId)				

			seqIdToAdd.append(maxMemberShipIds)

		seqsToAdd = list(itertools.product(*seqIdToAdd))
		for seqsAdded in seqsToAdd:
			finalClustersComplete.append(seqs  + list(seqsAdded))
	return(finalClustersComplete)

def sortClusters(clusters, x_coords):

	compactness ={}
	for cluster in clusters:
		val = 0.0
		for i in range(len(cluster)):
			for j in range(i, len(cluster)):
				val += compute_distance(x_coords[cluster[i]], x_coords[cluster[j]])
		compactness[val] = cluster

	sortedCluster = []
	for i in sorted (compactness.keys()): 
		#print (i, compactness[i])
		sortedCluster.append(compactness[i])

	return sortedCluster


def plotPoint(x_coords, X, finalClusterCentres, finalClusters, completeCluster):
	x_abs = []
	y_ord = []
	y_kmeans = []

	for x in X:
		x_abs.append(x_coords[x][0])
		y_ord.append(x_coords[x][1])
	#plt.scatter(x_abs, y_ord)

	x_center_abs = []
	y_center_ord = []	

	listCentre = []
	for center, coord in finalClusterCentres.items():
		x_center_abs.append(coord[0])
		y_center_ord.append(coord[1])
		listCentre.append(center)
	plt.scatter(x_center_abs, y_center_ord, c='black', s=200, alpha=0.5)

	x_abs_final = []
	y_ord_final = []
	points = []
	for x in X:
		for k, seq in finalClusters.items():
			if x in completeCluster:
				x_abs_final.append(x_coords[x][0])
				y_ord_final.append(x_coords[x][1])
				points.append([x_coords[x][0], x_coords[x][1]])
				y_kmeans.append(len(listCentre)+1)
				break
			if x in seq:
				y_kmeans.append(listCentre.index(k))

	points = np.asarray(points)
	hull = ConvexHull(points)
	for simplex in hull.simplices:
		plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

	#plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
	#plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
	
	plt.scatter(x_abs, y_ord, s=50, cmap='gray')
	#plt.scatter(x_abs, y_ord, c=y_kmeans, s=50, cmap='gray')
	plt.show()

def compute_alignment_muscle(nt_file, aa_file, nt_alignment, aa_alignment):

	command = "./muscle3.8.31_i86linux32 -in "  + aa_file + " -out " +  aa_alignment
	os.system(command)	

	nt_dict = {}
	for record in SeqIO.parse(nt_file, "fasta"):
		nt_dict[record.id] = str(record.seq)

	aa_aln_dict = {}
	for record in SeqIO.parse(aa_alignment, "fasta"):
		aa_aln_dict[record.id] = str(record.seq)		
	

	file = open(nt_alignment, "w")
	for id in aa_aln_dict.keys():
		back_translate_seq = ""
		aa_aln_seq = aa_aln_dict[id]
		nt_seq = nt_dict[id]
		nb_gap = 0
		for pos in range(len(aa_aln_seq)):
			if aa_aln_seq[pos] == "-":
				back_translate_seq += "---"
				nb_gap += 1
			else:
				pos_to_use = pos - nb_gap
				back_translate_seq += nt_seq[3*pos_to_use:3*pos_to_use+3]
		file.write(">" + id + "\n")
		file.write(back_translate_seq + "\n")

	file.close()



def createSeqForCluster(clusters, microalignmentFile,mappingGeneToTranscript, sub_path):
	try:		
		
		mappingGeneToTranscriptDetail = detailGeneTranscript(mappingGeneToTranscript)		

		path =  sys.path[0] + "/clusters_" + sub_path + "/fuzzyCMeans/"
		filename = microalignmentFile.split("/")[-1]
		filename = filename.split("_")[0]
		
		geneToSpecieFIle = open(sys.path[0] + "/geneToSpecie/cleanTree/" + filename + "_cleanoutput.out", "r")

		mapping_file = open("tmp/" + filename + "_mapping_species.txt", "w")
		mapping_transcript_file = open("tmp/" + filename + "_mapping_transcript.txt", "w")

		filename = filename + "_" + "cluster_"

		clusterFile = path + filename

		sequences = {}
		speciesOfGene = {}

		lines = geneToSpecieFIle.readlines()
		for line in lines:
			line = line.strip()
			tmp = line.split(" ")
			specieName = tmp[0]
			specieName = specieName.replace(" ", "")
			specieName = specieName.replace("/", "")
			specieName = specieName.replace("_", "")		
			specieName = specieName.replace("-", "")	
			specieName = specieName.upper()	
			mapping_file.write(specieName + "\t" + tmp[0] + "\n")
			speciesOfGene[tmp[1]] = specieName

		

		for record in SeqIO.parse(microalignmentFile, "fasta"):
			sequences[record.id] = record.seq

		i = 1
		filenames = [path]
			

		for seqId in mappingGeneToTranscriptDetail.keys():
			sequenceID = seqId + "" + mappingGeneToTranscriptDetail[seqId] + "_" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]]			
			mapping_transcript_file.write(sequenceID + "\t" + seqId + "\t" + mappingGeneToTranscriptDetail[seqId] +  "\t" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]] +"\n")
		mapping_file.close()
		mapping_transcript_file.close()
		#le programme enregistre uniquement les id des seq qui sont retenus par le master cluster
		for cluster in clusters:
			seqOfThisCLuster = []	
			AASeqOfThisCLuster = []						
			for seqId in cluster:
				sequence = str(sequences[seqId])
				sequence = Seq(sequence.replace("-", ""))
				#print(speciesOfGene)
				#print(mappingGeneToTranscriptDetail[seqId])
				#print(speciesOfGene[mappingGeneToTranscriptDetail[seqId]])
				#Ajouter for each cds of each gene
				sequenceID = seqId + "" + mappingGeneToTranscriptDetail[seqId] + "_" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]]
				
				#mapping_transcript_file.write(sequenceID + "\t" + seqId + "\t" + mappingGeneToTranscriptDetail[seqId] +  "\t" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]] +"\n")
				#print(sequenceID + "\t" + seqId + "\t" + mappingGeneToTranscriptDetail[seqId] +  "\t" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]])

				record = SeqRecord(sequence , sequenceID, '', '')
				seqOfThisCLuster.append(record)

				AARecord = SeqRecord(translate(sequence) , sequenceID, '', '')
				AASeqOfThisCLuster.append(AARecord)				
			
			filename_tmp = clusterFile + str(i) + ".fasta"		
			SeqIO.write(seqOfThisCLuster, filename_tmp, "fasta")

			filename_tmp2 = clusterFile + str(i) + "_seqAA.fasta"		
			SeqIO.write(AASeqOfThisCLuster, filename_tmp2, "fasta")

			nt_alignment = clusterFile + str(i) + "_NT.fasta"
			aa_alignment = clusterFile + str(i) + "_AA.fasta"

			compute_alignment_muscle(filename_tmp, filename_tmp2, nt_alignment, aa_alignment)

			"""
			command = "java -jar macse_v2.03.jar -prog alignSequences -seq "  + filename_tmp 
			os.system(command)	
			"""

			filenames.append(filename_tmp)

			sequences_nt = []
			for record in SeqIO.parse(nt_alignment, "fasta"):
				sequence = str(record.seq)
				sequence = Seq(sequence.replace("!", "A"))			
				sequences_nt.append(SeqRecord(sequence, record.id, '', ''))  

			SeqIO.write(sequences_nt, nt_alignment, "fasta")

			"""
			command2 = "./treebest best " + nt_alignment + " > " +  clusterFile + str(i) + "_.nhx -f ressources/ensemblSpecies.tree"
			print (command2)
			os.system(command2)

			command3 = "./treebest root " + clusterFile + str(i) + "_.nhx" + " > " +  clusterFile + str(i) + "_root.nhx"
			print (command3)
			os.system(command3)
			"""
			command2 = "./treebest best " + nt_alignment + " > " +  clusterFile + str(i) + "_root.nhx -f ressources/ensemblSpecies.tree"
			#print (command2)
			os.system(command2)

			#command3 = "./treebest root " + clusterFile + str(i) + "_.nhx" + " > " +  clusterFile + str(i) + "_root.nhx"
			#print (command3)
			#os.system(command3)

			i +=1
		

		return filenames
	except Exception, e:
		print(e)
		#mapping_file.close()
		#mapping_transcript_file.close()
		return []

def writeClusters(finalClusters, finalClustersclustersComplete, microalignmentFile, sub_path):
	path =  sys.path[0] + "/clusters_" + sub_path +"/clusters/"
	filename = microalignmentFile.split("/")[-1]

	filename = filename.split("_")[0]

	subClusters = open(path + filename + "_subClusters.txt", "w")

	completeClusters = open(path + filename + "_completeClusters.txt", "w")

	i = 0
	for k, v in finalClusters.items():
		subClusters.write(">cluster_" + str(i) + "\n")
		subClusters.write("\t".join(v) + "\n")
		i +=1
		

	j = 0
	for cluster in finalClustersclustersComplete:
		completeClusters.write(">completeCluster_" + str(j) + "\n")
		completeClusters.write("\t".join(cluster) + "\n")
		j +=1
	subClusters.close()
	completeClusters.close()
# A Simple Merge based O(n) Python 3 solution  
# to find median of two sorted lists 
  
# This function returns median of ar1[] and ar2[]. 
# Assumptions in this function: 
# Both ar1[] and ar2[] are sorted arrays 
# Both have n elements 
def getMedian( ar1, ar2): 
	result = []
	for i in range(len(ar1)):
		result.append((ar1[i] + ar2[i])/2.0)
	return result

def upgma(x_coords, finalClustersclustersComplete, microalignmentFile, sub_path):
	#print(x_coords, finalClustersclustersComplete)
	path =  "clusters_" + sub_path + "/fuzzyCMeans/"
	
	filename = microalignmentFile.split("/")[-1]

	filename = filename.split("_")[0]
	#print(finalClustersclustersComplete, "1111")
	tree = ""
	iteration = 0
	for cluster in finalClustersclustersComplete:
		iteration += 1

		file = open(path + filename + "_upgma_" + str(iteration)  + "_root.nhx", "w")
		seqIds = cluster
		while(len(seqIds)>1):
			min = 1.1
			seq_i = ""
			seq_j = ""
			nb_elt = len(seqIds)

			for i in range(nb_elt):
				for j in range(i):
					seq1 = seqIds[i]
					seq2 = seqIds[j]
					if compute_distance(x_coords[seq1], x_coords[seq2]) < min and seq1 in seqIds and seq2 in seqIds:
						min = compute_distance(x_coords[seq1], x_coords[seq2])
						seq_i = seq1
						seq_j = seq2 
			seqIds.remove(seq_i)
			seqIds.remove(seq_j)
			seqIds.append("(" + seq_i + "," + seq_j + ")")
			x_coords["(" + seq_i + "," + seq_j + ")"] =   getMedian(x_coords[seq1], x_coords[seq2])
		
		tree = seqIds[0] + ";"
		file.write(tree)
		file.close()

def upgma_2(x_coords, finalClustersclustersComplete, microalignmentFile):
	#print(x_coords, finalClustersclustersComplete)
	path =  "clusters_50_50/fuzzyCMeans/"
	
	filename = microalignmentFile.split("/")[-1]

	filename = filename.split("_")[0]

	tree = ""
	iteration = 0
	for cluster in finalClustersclustersComplete:
		iteration += 1

		file = open(path + filename + "_transcrit_" + str(iteration)  + "_root.nhx", "w")
		seqIds = cluster
		while(len(seqIds)>1):
			min = 1.1
			seq_i = ""
			seq_j = ""
			nb_elt = len(seqIds)
			for i in range(nb_elt):
				for j in range(i):
					seq1 = seqIds[i]
					seq2 = seqIds[j]

					if compute_distance(x_coords[seq1], x_coords[seq2]) < min and seq1 in seqIds and seq2 in seqIds:
						#print("___________", seq1, seq2)
						min = compute_distance(x_coords[seq1], x_coords[seq2])
						seq_i = seq1
						seq_j = seq2 

			seqIds.remove(seq_i)
			seqIds.remove(seq_j)
			seqIds.append("(" + seq_i + "," + seq_j + ")")
			x_coords["(" + seq_i + "," + seq_j + ")"] =   getMedian(x_coords[seq1], x_coords[seq2])
		
		tree = seqIds[0] + ";"
		t = Tree(tree)

		file.write(tree)
		file.close()


def NJ(x_coords, finalClustersclustersComplete, microalignmentFile, sub_path):
	#print(x_coords, finalClustersclustersComplete)
	path =  "clusters_" + sub_path + "/fuzzyCMeans/"
	
	filename = microalignmentFile.split("/")[-1]

	filename = filename.split("_")[0]

	tree = ""
	iteration = 0
	#print(finalClustersclustersComplete, "22222")
	for cluster in finalClustersclustersComplete:
		iteration += 1

		#file = open(path + filename + "_nj_" + str(iteration)  + "_.nhx", "w")
		file = open(path + filename + "_nj_" + str(iteration)  + "_root.nhx", "w")
		seqIds = cluster
		distanceMatrix = []
		header = []
		for seq1 in seqIds:
			line = []
			for seq2 in seqIds:
				line.append(compute_distance(x_coords[seq1], x_coords[seq2]))
			distanceMatrix.append(line)
			header.append(seq1)

		dm = DistanceMatrix(distanceMatrix, header)
		tree = nj(dm)

		newick_str = nj(dm, result_constructor=str)
		file.write(newick_str)
		file.close()
		"""
		command = "./treebest root " + path + filename + "_nj_" + str(iteration)  + "_.nhx" + " > " +  path + filename + "_nj_" + str(iteration)  + "_root.nhx"
		print (command)
		os.system(command)		
		"""

def main_fuzzy_cMeans(X, x_coords, mappingGeneToTranscript, microalignmentFile, sub_path):	

	if len(X) <= 2:
		filenames = createSeqForCluster([X], microalignmentFile, mappingGeneToTranscript, sub_path)
		return (filenames)
	#print(X, x_coords)
	fuzzifier = 2
	epsilon = 0.0001


	#x_coords ={'x10': [5, 6], 'x8': [6, 10], 'x9': [2, 6], 'x11': [2, 2], 'x2': [7, 8], 'x3': [9, 9], 'x12': [10, 9], 'x1': [9, 4], 'x6': [8, 9], 'x7': [8, 6], 'x4': [2, 10], 'x5': [5, 6]}
	
	n = len(X)
	#print(X)

	Cmin, Cmax = chooseCminCmax(n, mappingGeneToTranscript)	

	min_value_of_validity_index = 99999.0
	optimalCLuster = {}
	allCLuster = {}
	sep_c_max = 0.0

	clusters = {}
	clusterCentersIndex = []
	for i in range(Cmin):
		clusters["v"+str(i)] = ""
		clusterCentersIndex.append("v"+str(i))

	clusterCenters = initializeClusterCenters(clusters, X, x_coords)		


	for c in range(Cmin, Cmax+1):
		#print (c)
		#print (x_coords)		
		#clusters = {}
		#clusterCentersIndex = []
		#for i in range(c):
		#	clusters["v"+str(i)] = ""
		#	clusterCentersIndex.append("v"+str(i))

		#clusterCenters = initializeClusterCenters(clusters, X, x_coords)		

		#print ("=======================",clusterCenters)
		U = np.zeros(shape=(len(X),len(clusterCentersIndex)))


		clusterCenters, menbershipMatrix = FCM(U, fuzzifier, X, x_coords, clusterCentersIndex, clusterCenters, epsilon)		
		#print (menbershipMatrix)
		#print(clusterCenters, menbershipMatrix)
		

		clusters = extractCluster(menbershipMatrix, clusterCenters, clusterCentersIndex, X, mappingGeneToTranscript)	

		##print (clusterCenters)
		#print(clusters)	

		#print (U)	
		scat_c, sep_c = validityIndex(X, x_coords, clusterCentersIndex, clusterCenters, U)
		allCLuster[c] = [copy.deepcopy(scat_c), copy.deepcopy(sep_c), copy.deepcopy(clusters), copy.deepcopy(clusterCenters), copy.deepcopy(clusterCentersIndex), copy.deepcopy(U)]
		#sep_c_max = sep_c
		if sep_c > sep_c_max:
			sep_c_max = sep_c
		#print ("===",allCLuster, "===")
		#print (c,scat_c, sep_c)
		

		worstClusterToSplit = findWorstClusterToSplit(U, clusters, X, x_coords, clusterCentersIndex, clusterCenters)
		if len(worstClusterToSplit) == 0:
			break
		else:
			clusterCenterToSplit = worstClusterToSplit[0]
			sequencesToSplit = worstClusterToSplit[1]

		#print(clusterCenterToSplit, sequencesToSplit)		
				
		v_i1, coord_v_i1, v_i2, coord_v_i2, E1, E2  = findNextTwoCenters(clusterCenterToSplit, sequencesToSplit, clusterCentersIndex, clusterCenters, U, x_coords)

		#print(v_i1, coord_v_i1, v_i2, coord_v_i2)
		#exit()	
		clusterCenters[v_i1] = coord_v_i1
		clusterCenters[v_i2] = coord_v_i2
		#x_coords[v_i1] = coord_v_i1
		#x_coords[v_i2] = coord_v_i2
		clusterCentersIndex.append(v_i1)
		clusterCentersIndex.append(v_i2)

		clusterCentersIndex.remove(clusterCenterToSplit)
		del clusterCenters[clusterCenterToSplit]
		#print (clusters)		
		del clusters[clusterCenterToSplit]
		clusters[v_i1] = E1
		clusters[v_i2] = E2

		#print("\n___________________")
	

	for c in range(Cmin, Cmax+1):	
		#print("\n=============================")
		scat_c = allCLuster[c][0]
		sep_c = allCLuster[c][1]

		validity_index = scat_c + sep_c/sep_c_max
		#print(validity_index, allCLuster[c][2], allCLuster[c][3])
		
		if validity_index < min_value_of_validity_index:
			min_value_of_validity_index = validity_index
			optimalCLuster = allCLuster[c]

	#print("\n________________________________")
	#print(min_value_of_validity_index)
	#print(optimalCLuster)
	finalClusters = optimalCLuster[2]
	finalClusterCentres = optimalCLuster[3]
	finalclusterCentersIndex = optimalCLuster[4]
	U = optimalCLuster[5]
	#print(finalClusters)
	#print(finalClusterCentres)
	#print(U)
	finalClustersclustersComplete = selectAndCompleteTranscriptTree(finalClusters, finalClusterCentres, finalclusterCentersIndex, mappingGeneToTranscript, X, x_coords, U)
	
	#print("\n_____________________________________________")
	#print(finalClusters)
	#print(finalClustersclustersComplete)
	#"""
	writeClusters(finalClusters, finalClustersclustersComplete, microalignmentFile, sub_path)
	
	"""
	to have an overview of visual representation uncomment this two lines
	for completeCluster in finalClustersclustersComplete:
		plotPoint(x_coords, X, finalClusterCentres, finalClusters, completeCluster)
	"""
	compactness = sortClusters(finalClustersclustersComplete, x_coords)
	clustersList = []
	for k, v in finalClusters.items():
		clustersList.append(v)
	upgma_2(x_coords, clustersList, microalignmentFile)

	filenames = createSeqForCluster(compactness, microalignmentFile, mappingGeneToTranscript, sub_path)
	#print(x_coords, finalClustersclustersComplete, microalignmentFile)
	#upgma(x_coords, copy.deepcopy(finalClustersclustersComplete), microalignmentFile)
	
	#print(x_coords, finalClustersclustersComplete, microalignmentFile)
	NJ(x_coords, finalClustersclustersComplete, microalignmentFile, sub_path)

	return filenames
	#print("\n_____________________________________________")

if __name__ == "__main__":
	#upgma({'ENSCABT00000028818':[-0.87065561], 'ENSCABT00000028809': [-0.89092845], 'ENSCABT00000028814': [-0.91891197], 'ENSGAGT00000004415': [-0.80850891], 'ENSCABT00000028813': [-0.91891197], 'ENSGAGT00000004427': [-0.80038165]}, [['ENSCABT00000028818', 'ENSGAGT00000004415'], ['ENSGAGT00000004427', 'ENSCABT00000028814']])

	#exit()
	X = ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12"]
	mappingGeneToTranscript = {"g1":["x1", "x2"], "g2":["x3"], "g3":["x4", "x5"],  "g4":["x6", "x7"],  "g5":["x8", "x9", "x10", "x11"],  "g6":["x12"]}
	x_coords ={'x10': [5, 6], 'x8': [6, 10], 'x9': [2, 6], 'x11': [2, 2], 'x2': [7, 8], 'x3': [9, 9], 'x12': [10, 9], 'x1': [9, 4], 'x6': [8, 9], 'x7': [8, 6], 'x4': [2, 10], 'x5': [5, 6]}

	x_coords = {}
	for  xi in X:
		x_coords[xi] = [random.randint(1,10), random.randint(1,10)]

	main_fuzzy_cMeans(X, x_coords, mappingGeneToTranscript, "")
