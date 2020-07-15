#!/usr/bin/python

'''
Author : Marie Degen

Labeled all the exons with his corresponding alternative splicing event from a gene family. 
Example : python alternativeSplicing.py spliceGraph/FAM86_spliceGraphe.txt macroalignment/FAM86_macroalignment.txt microalignment/FAM86_microalignment.fasta

'''

import sys
import argparse
import requests
import os
from os import getcwd
from collections import Counter
import operator
from Bio import SeqIO
import time

def build_arg_parser(path):
	parser = argparse.ArgumentParser(description="Ensembl gene tree parsor program parameters")   
	parser.add_argument('-s', '--save_path', default =  path + '/labelExonsMaker/')
	parser.add_argument('-mi', '--microalignment',  required=True)
	parser.add_argument('-ma', '--macroalignment',  required=True)
	parser.add_argument('-i', '--spliceGraph', required=True)
	return parser


def  microalignmentIntoDictionnary(file, dictSequence):
	#microalignment into dictionary
	for record in SeqIO.parse(file, "fasta"):
		dictSequence[record.id] = str(record.seq)

		
def macroalignmentIntoList(file, macro):		
	#macroalignment into dictionary
	tmp = file.split("/")
	tmp = tmp[-1]
	tmp = tmp.split("_")[0]
	source2target = open("initialSource/" +  tmp + "_initialsource2target.txt", "r")
	lines = source2target.readlines()
	transcripts = []
	genes =[]
	for line in lines:
		line = line.replace("\n", "")
		parts = line.split(" ")
		transcripts.append(parts[0])
		genes.append(parts[1])

	exonList = []
	with open(file,'r') as macro_file:
		for macro_line in macro_file: 
			#exon ID 
			if macro_line.startswith(">"):
				if not exonList :
					continue
				else:
					macro.append(exonList)
					exonList = []				
			else: 
				#get the id 
				s = macro_line.split(":")
				if len(s)>=2:
					i = s[0]
					l = s[1]
					#get the type
					if i in genes:
						continue
					elif i in transcripts:
						start = l.split("-")[0]
						end = l.split("-")[1]
						end = end.replace("\n", "")
						end = int(end)
						if end != 0:
							end -=1
						exonList.append([i, start, end])

		macro.append(exonList)
		exonList = []	
	

def getPhaseLength(seqTranscrit, res):
	#get the phase of the exon 
	newSeqTranscrit = seqTranscrit.replace("-", "")
	
	phase = abs(len(seqTranscrit) - len(newSeqTranscrit))%3
	
	#add the length
	res.append(len(newSeqTranscrit))
	#res.append(len(seqTranscrit))
	
	if phase == 0:
		res.append(0)
	elif phase == 1:
		res.append(1)
	elif phase == 2:
		res.append(2)
	else: 
		print ("Error in the phase calcul")
	
	return res

def getFinalSequence(idRefExon,startRef, start, endRef, end, dictSequence):
	sequenceRefExon = dictSequence.get(idRefExon)
	counterRef = 0 
	finalSequenceRef = ""
	
	maxmin = [startRef, start,endRef, end]
	maximum = max(maxmin)
	minimum = min(maxmin)
	#retrieve refSequence	
	#time.sleep(1)
	for q in sequenceRefExon:
		if counterRef < maximum: 
			if q == '-' and counterRef < minimum :
				pass 
			elif q != '-' and counterRef >= minimum: 
				counterRef += 1
				finalSequenceRef += q
			elif q == '-' and counterRef >= minimum:
				counterRef += 1
				finalSequenceRef += q
			else :
				counterRef += 1
	#print(idRefExon, finalSequenceRef)
	#print("-", idRefExon, startRef, start,endRef, end, finalSequenceRef) 
	return finalSequenceRef


def getFinalSequence_old(idRefExon,startRef, start, endRef, end, dictSequence):
	sequenceRefExon = dictSequence.get(idRefExon)
	counterRef = 0 
	finalSequenceRef = ""
	
	maxmin = [startRef, start,endRef, end]
	maximum = max(maxmin)
	minimum = min(maxmin)
	#retrieve refSequence	
	time.sleep(1)
	for q in sequenceRefExon:
		if counterRef < maximum: 
			if q == '-' and counterRef < minimum :
				pass 
			elif q != '-' and counterRef >= minimum: 
				counterRef += 1
				finalSequenceRef += q
			elif q == '-' and counterRef >= minimum:
				counterRef += 1
				finalSequenceRef += q
			else :
				counterRef += 1
	#print(idRefExon, finalSequenceRef)
	#print("-", idRefExon, startRef, start,endRef, end, finalSequenceRef) 
	return finalSequenceRef



def defineEvent35prime(idRefExon, idExonMacro, startRef, start, endRef, end, constitutive, exonskipping, dictSequence):
	finalSequenceRef = ""
	finalSequenceExon = ""
	tprimeextension = 0 
	tprimeretention = 0 
	cprimeextension = 0 
	cprimeretention = 0 
	nosplicingevent = 0
	res = [] 
	
	tmp =  dictSequence.get(idRefExon)
	tmp = tmp[startRef:endRef]
	gapsLen = tmp.count("-")

	finalSequenceRef = getFinalSequence(idRefExon, startRef, start, endRef+gapsLen, end, dictSequence)
	#print (finalSequenceRef)
	#exit()
				
	finalSequenceExon = getFinalSequence(idExonMacro, startRef, start, endRef+gapsLen, end, dictSequence)

	if constitutive == True: 
		res.append('preserved')
	elif exonskipping == True: 
		res.append('exon skipping')
	else : 
		#5'extension
		for o,p in zip(finalSequenceRef, finalSequenceExon):
			if o == '-' and p == '-':
				continue
			elif o != '-' and p == '-':
				cprimeretention += 1
				break
			elif p != '-' and o == '-':		
				cprimeextension += 1
				break
			else :
				break

		#3'extension
		for o,p in zip(reversed(finalSequenceRef), reversed(finalSequenceExon)):
			if o == '-' and p == '-':
				continue
			elif o != '-' and p == '-':
				tprimeretention += 1
				break
			elif p != '-' and o == '-':		
				tprimeextension += 1
				break
			else:
				break

		#print str(cprimeretention), str(cprimeextension),str(tprimeretention),str(tprimeextension)
		
		if cprimeretention != 0  and tprimeretention != 0 :
			res.append('3 retention and 5 retention')
		elif cprimeextension != 0  and tprimeextension != 0 :
			res.append('3 extension and 5 extension')
		elif cprimeretention != 0  and tprimeextension != 0 :
			res.append('3 extension and 5 retention')
		elif cprimeextension != 0  and tprimeretention != 0 :
			res.append('3 retention and 5 extension')	
		elif cprimeretention != 0 : 
			res.append('5 retention')
		elif cprimeextension != 0 : 
			res.append('5 extension')
		elif tprimeretention != 0 : 
			res.append('3 retention')
		elif tprimeextension != 0 : 
			res.append('3 extension')
		else: 
			res.append('no splicing event')
		
		
	#get the phase and length of the exon 	
	res =  getPhaseLength(finalSequenceExon, res)
	return res


def main_alternativeSplicing(save_path, macroalignment, microalignment, splicegraph):
	refID = ""
	refPosition = ""
	exonID = ""
	dictSequence = {}
	macro = []
	refExonList = []
	microalignmentIntoDictionnary(microalignment, dictSequence)
	macroalignmentIntoList(macroalignment, macro)

	#splicegraph into a list
	with open(splicegraph) as splice_file:
		namefile = splicegraph.split("/")[-1]
		namefile = namefile.split(".")[0]
		namefile = namefile + "_labeledExons.out"
		completeName = os.path.join(save_path, namefile)         
		outputfile = open(completeName, "w")
		indexMacro = 0 
		for line in splice_file :		
			#print (line)
			#exon ID 
			if line.startswith(">"):
				continue
			#Transcrit ID and position 
			else:
				if len(line.split(":"))>=2:
					refID = line.split(":")[0]
					refPosition = line.split(":")[1]
					startRef = int(refPosition.split("-")[0])
					endRef = refPosition.split("-")[1]
					endRef = endRef.replace("\n", "")
					endRef = int(endRef)
					refExonList.append([refID, startRef, endRef])
					
	#compare between splicegraph and macro/micro
	#print(macro)
	for n in range(len(refExonList)):
		outputfile.write(">" + str(n) + "\n")
		#loop on macro alignment
		for m in range(len(macro[n])):
			exonskipping = False
			constitutive = False 
			#conserved exon
			if int(refExonList[n][1]) == int(macro[n][m][1]) and int(refExonList[n][2]) == int(macro[n][m][2]):
				constitutive = True
				#print("constitutive")
				res = defineEvent35prime(refExonList[n][0],macro[n][m][0],int(refExonList[n][1]), int(macro[n][m][1]), int(refExonList[n][2]), int(macro[n][m][2]), constitutive, exonskipping, dictSequence)
				outputfile.write(macro[n][m][0] + ":"+ str(res) + "\n")
			#exon skipping = cassette exon 
			elif int(macro[n][m][1]) == 0 and int(macro[n][m][2]) == 0:
				exonskipping = True
				#print("exonskipping")
				res = defineEvent35prime(refExonList[n][0],macro[n][m][0],int(refExonList[n][1]), int(macro[n][m][1]), int(refExonList[n][2]), int(macro[n][m][2]), constitutive, exonskipping, dictSequence)
				outputfile.write(macro[n][m][0] + ":"+ str(res) + "\n")
			#3' extension retention or 5' extension retention
			elif int(refExonList[n][1]) != int(macro[n][m][1]) or int(refExonList[n][2]) != int(macro[n][m][2]):
				#check into microalignment call function 
				#print("5prime, 3'")
				res = defineEvent35prime(refExonList[n][0],macro[n][m][0],int(refExonList[n][1]), int(macro[n][m][1]), int(refExonList[n][2]), int(macro[n][m][2]), constitutive, exonskipping, dictSequence)
				outputfile.write(macro[n][m][0] + ":"+ str(res) + "\n")
			else:
				outputfile.write(macro[n][m][0] + ":another event \n")
		#exit()
	outputfile.close()
	return completeName


def alternativeSplicing_aux(macroalignment, microalignment, splicegraph):
	path =  sys.path[0]      
	save_path  = path + "/labelExonsMaker/"	
	return main_alternativeSplicing(save_path, macroalignment, microalignment, splicegraph)


if __name__ == "__main__":    
	path =  sys.path[0]
	parser = build_arg_parser(path)
	arg = parser.parse_args()        
	save_path  = arg.save_path 
	macroalignment = arg.macroalignment
	microalignment = arg.microalignment
	splicegraph = arg.spliceGraph
	main_alternativeSplicing(save_path, macroalignment, microalignment, splicegraph)