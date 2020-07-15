#!/usr/bin/python

'''
Author : Marie Degen


Cluster the transcrits from a genes family
Example : python retrieveSeqClusters.py clusters/FAM86_clusters.out microaligment/FAM86_microalignment.fasta

'''

import sys
import argparse
import os
from Bio import SeqIO
import pickle
import shutil


def build_arg_parser(path):
	parser = argparse.ArgumentParser(description="Label exons")   
	parser.add_argument('-s', '--save_path', default =  path + '/clusters/sequencesClusters/')
	parser.add_argument('-c', '--cluster_file',  required=True)
	parser.add_argument('-m', '--microalignment', required=True )
	return parser


def  microalignmentIntoDictionnary(file, dictSequence):
	#microalignment into dictionary
	for record in SeqIO.parse(file, "fasta"):
		dictSequence[record.id] = str(record.seq)
		

def getCluster(file, dictSequence, save_path):
	with open(file,'r') as cluster_file:
		for cluster_line in cluster_file:
			cluster_line = cluster_line.replace("[","")
			cluster_line = cluster_line.replace("]","")
			cluster_line = cluster_line.replace("'","")
			
			 
			nameDirectory = file.split("/")[-1].split(".")[0]
			completenameDire = os.path.join(save_path, nameDirectory)   
			
			namephylip = "phylipFormat"
			completenamePhyl = os.path.join(completenameDire, namephylip) 
			
			if not os.path.exists(completenameDire):
				os.mkdir(completenameDire)
			else:
			    shutil.rmtree(completenameDire)           
			    os.mkdir(completenameDire)

			if not os.path.exists(completenamePhyl):
				os.mkdir(completenamePhyl)
			else:
			    shutil.rmtree(completenameDire)           
			    os.mkdir(completenameDire)			
			i = 0
			cluster = cluster_line.split(",")
			while i < len(cluster):
				namefile = "group" + str(i) + ".fasta"
				completeName = os.path.join(completenameDire, namefile)         
				outputfile = open(completeName, "w")
				
				cluster_line = cluster[i].replace("(","")
				cluster_line =  cluster_line.replace(")","")
				
				couple = cluster_line.split("-")
				j = 0 
				while j < len(couple):
					transcriptID = couple[j].split("_")[0]
					transcriptID = transcriptID.replace(" ","")
					outputfile.write(">" + transcriptID + "\n")
					outputfile.write(dictSequence.get(transcriptID) + "\n")
					j+=1
					
				i+=1
	return completenameDire


def main_retrieve_cluster(save_path, cluster_file, microalignment):
	dictSequence = {}
	groups =[]
	microalignmentIntoDictionnary(microalignment, dictSequence)
	return getCluster(cluster_file, dictSequence, save_path)


def main_retrive_cluster_aux(cluster_file, microalignment):
	path =  sys.path[0]      
	save_path  = path + "/clusters/sequencesClusters/"
	return main_retrieve_cluster(save_path, cluster_file, microalignment)


if __name__ == "__main__":    
	path =  sys.path[0]
	parser = build_arg_parser(path)
	arg = parser.parse_args()        
	cluster_file  = arg.cluster_file 
	microalignment = arg.microalignment
	save_path = arg.save_path
	main_retrieve_cluster(save_path, cluster_file, microalignment)
