#!/usr/bin/python

'''
Author : Marie Degen


Score all the transcrits of a family genes
Example : python -W ignore matrixScore.py labelExons/FAM86_labeledExons.out initialSource/FAM86_rcg17_initialsource2target.txt

'''

import sys
import argparse
import os
from os import getcwd
from collections import Counter
import operator
import numpy as np
import pandas as pd


def build_arg_parser(path):
	parser = argparse.ArgumentParser(description="Label exons")   
	parser.add_argument('-s', '--save_path', default =  path + '/similarityScores/')
	parser.add_argument('-l', '--labeledExons',  required=True)
	parser.add_argument('-i', '--initialsource2target',  required=True)
	parser.add_argument('-t', '--type',  default = "microalignment")

	return parser

def scoreLength(length1, length2):
	if int(length1) == 0 and int(length2) == 0:
		return 1
	else: 
		return 1-(abs(int(length1) - int(length2)) / (int(length1) + int(length2)))

def scorePhase(phase1, phase2):
	return 1 if int(phase1) == int(phase2) else 0

def scoreEvent(transcrit1, transcrit2):
	listeventtranscrit1 = []
	listeventtranscrit2 = []
	position1 = []
	position2 = []
	nbtotalevent = 0
	nbeventcommon = 0 
	exon = 1
	exonskipping = False
	
	for exon1 in transcrit1:
		if not exon1.isdigit():
			listeventtranscrit1.append(exon1)
		else :
			position1.append(exon1)
			
	for exon2 in transcrit2:
		if not exon2.isdigit():
			listeventtranscrit2.append(exon2)
		else :
			position2.append(exon2)

	for i in listeventtranscrit1:
		for j in listeventtranscrit2:
			if i == "'exonskipping'" and j == "'exonskipping'":
				exon = 0
				exonskipping = True
				return exon,exonskipping, 0.0
			elif i == "'exonskipping'" or j == "'exonskipping'":
				exon = 1
				exonskipping = True
				return exon,exonskipping, 0.0
			elif i == j : 
				if i  == "'internalgap<30'" and j  == "'internalgap<30'" or i == "'internal gap >= 30'" and j == "'internal gap >= 30'" : 
					if position1[0] == position2[0]:
						nbeventcommon += 1
						nbtotalevent += 1
					else:
						nbtotalevent += 1
				else: 
					nbeventcommon += 1
					nbtotalevent += 1
			else: 
				nbtotalevent += 1	
				
	return exon,exonskipping, float(nbeventcommon / nbtotalevent)	

def separateEventLengthPhase(exon):
	l = len(exon)-2
	lengthPhase = exon[-2:]
	event = exon[:l]
	return event, lengthPhase
		
def formatTranscrit(ID, transcrit):
	transcrit = transcrit.replace("\n","")
	transcrit = transcrit.replace("["," ")
	transcrit = transcrit.replace("]","")
	if len(transcrit.split(",")) == 3:
		return [ID,transcrit.split(",")[0].replace(" ",""), transcrit.split(",")[1].replace(" ",""),transcrit.split(",")[2].replace(" ","")]
	elif len(transcrit.split(",")) == 4:
		return [ID,transcrit.split(",")[0].replace(" ",""), transcrit.split(",")[1].replace(" ",""),transcrit.split(",")[2].replace(" ",""), transcrit.split(",")[3].replace(" ","")]
	elif len(transcrit.split(",")) == 5:
		return [ID,transcrit.split(",")[0].replace(" ",""), transcrit.split(",")[1].replace(" ",""),transcrit.split(",")[2].replace(" ",""), transcrit.split(",")[3].replace(" ",""), transcrit.split(",")[4].replace(" ","")]
	elif len(transcrit.split(",")) == 6:
		return [ID,transcrit.split(",")[0].replace(" ",""), transcrit.split(",")[1].replace(" ",""),transcrit.split(",")[2].replace(" ",""), transcrit.split(",")[3].replace(" ",""), transcrit.split(",")[4].replace(" ",""), transcrit.split(",")[5].replace(" ","")]
	else:
		return [ID,transcrit.split(",")[0].replace(" ",""), transcrit.split(",")[1].replace(" ",""),transcrit.split(",")[2].replace(" ",""), transcrit.split(",")[3].replace(" ",""), transcrit.split(",")[4].replace(" ",""), transcrit.split(",")[5].replace(" ",""), transcrit.split(",")[6].replace(" ","")]

	
def compareTranscrit(transcrit1, transcrit2, alpha = 1, beta = 0, gamma = 1):
	scoreTotal = alpha +  beta + gamma
	sumexon = 0.0
	n = 1
	nbexontotal = 0 
	#print(transcrit1, transcrit2)
	for o in range(len(transcrit1)):
		score_event = 0
		score_length = 0 
		score_phase = 0
		if transcrit1[0] == transcrit2[0]:
			return 1
		else: 
			if int(n) == int(len(transcrit1)):
				break
			else:
				event1, lengthPhase1 = separateEventLengthPhase(transcrit1[n])
				event2, lengthPhase2 = separateEventLengthPhase(transcrit2[n])
				exon, exonskipping, score_event = scoreEvent(event1, event2)
				
				if not exonskipping : 
					score_length = scoreLength(lengthPhase1[0], lengthPhase2[0])
					score_phase = scorePhase(lengthPhase1[1], lengthPhase2[1])
					#print(lengthPhase1[1], lengthPhase2[1], score_phase)
					
				nbexontotal += exon

				sumexon += float(float(score_event*alpha) + float(score_length*beta) + float(score_phase*gamma)) / float(scoreTotal)
		n += 1
	
	
	return sumexon / nbexontotal


def main_matrixScore(save_path, labeledExons, initialsource2target, type):
	macro = []
	nbexon = 0
	exonList = []

	with open(labeledExons,'r') as label_file:
		for label_line in label_file: 
			#exon ID 
			if label_line.startswith(">"):
				nbtranscrit = 0
				nbexon += 1
				if not exonList :
					continue
				else:
					macro.append(exonList)
					exonList = []
			else: 
				nbtranscrit += 1
				#get the id 
				s = label_line.split(":")
				if len(s)>=2:
					exonList.append(formatTranscrit(s[0], s[1]))	
		macro.append(exonList)
		exonList = []
	transcritList = []
	#print macro
		

	for j in range(nbtranscrit):
		transcrit = []
		transcrit.append(macro[0][j][0])
		#print macro 
		for i in range(nbexon):
			temp = []
			for z in range(len(macro[i][j])):
				nb = z + 1
				if nb == len(macro[i][j]):
					break
				else:
					temp.append(macro[i][j][nb])
			transcrit.append(temp)
		transcritList.append(transcrit)


	#creation of the score matrix 	
	matriceScore = np.eye(nbtranscrit,nbtranscrit)
	names = []
	for k in range(len(transcritList)):
		names.append(transcritList[k][0])
		for l in range(len(transcritList)):
			res = compareTranscrit(transcritList[k],transcritList[l])
			matriceScore[k][l] = float("{0:.2f}".format(res))
			#matriceScore[l][k] = float("{0:.2f}".format(res))

	#get id transcrit with id gene
	with open(initialsource2target,'r') as target_file:
		o = 0 		
		for target_line in target_file: 
			if target_line.split(" ")[0] == names[o]:
				names[o] = target_line.split(" ")[0] + "_" + target_line.split(" ")[1]
				o += 1
	
	namefile = labeledExons.split("/")[-1]
	namefile = namefile.split("_")[0]
	namefile = type + "_" +namefile + "_score.csv"
	completeName = os.path.join(save_path, namefile)         

	df = pd.DataFrame(matriceScore, index=names, columns=names)

	with open(completeName, 'w') as outputfile:
		df.to_csv(outputfile)	
	return completeName

def main_matrixScore_aux(labeledExons, initialsource2target, type= "microalignment"):
	path =  sys.path[0]      
	save_path  = path + "/similarityScores/"
	return main_matrixScore(save_path, labeledExons, initialsource2target, type)

if __name__ == "__main__":    
	path =  sys.path[0]
	parser = build_arg_parser(path)
	arg = parser.parse_args()        
	save_path  = arg.save_path 
	initialsource2target = arg.initialsource2target
	labeledExons = arg.labeledExons
	type = arg.type
	main_matrixScore(save_path, labeledExons, initialsource2target, type)
