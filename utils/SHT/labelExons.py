#!/usr/bin/python

'''
Author : Marie Degen


Labeled the exons from the microalignment without the spliceGraphe, only base on the sequence
Example : python labelExons.py macroaligment/FAM86_macroalignment.txt microaligment/FAM86_microalignment.fasta


'''

import sys
import argparse
import os
from Bio import SeqIO



def build_arg_parser(path):
	parser = argparse.ArgumentParser(description="Label exons")   
	parser.add_argument('-s', '--save_path', default =  path + '/labelExons/')
	parser.add_argument('-mi', '--microalignment',  required=True)
	parser.add_argument('-ma', '--macroalignment',  required=True)
	return parser

'''
Function to put the microalignment into a dictionary, with the ID of the transcript as the key 
'''
def  microalignmentIntoDictionnary(file, dictSequence):
	#microalignment into dictionary
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


	transcriptID = []
	for record in SeqIO.parse(file, "fasta"):
		id = record.id
		if id in transcripts:
			dictSequence[record.id] = str(record.seq)
			transcriptID.append(record.id)
	"""
		id = record.id
		if id in genes:
			continue
		else:
			dictSequence[record.id] = str(record.seq)
			transcriptID.append(record.id)
	if len(transcriptID)>0:
		for i in range(len(dictSequence[transcriptID[0]])):
			flag = True
			for transcrit_id in transcriptID:
				if dictSequence[transcrit_id][i] != "-":
					flag = False
			if flag:
				for transcrit_id in transcriptID:
					seq = list(dictSequence[transcrit_id])
					seq[i] = "+"
					seq = "".join(seq)
					dictSequence[transcrit_id] = seq

	fileUpdate = open(file, "w") 					

	for transcrit_id in transcriptID:
		seq = dictSequence[transcrit_id]
		seq = seq.replace("+", "")		
		dictSequence[transcrit_id] = seq
		fileUpdate.write(">"+ transcrit_id+"\n")
		fileUpdate.write(str(seq) + "\n")
	"""
	


'''
Function to put the macroalignment into a array, id - start position - end position are added to the array 
'''		
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
'''
Return the maximum length of an a given exon
'''	
def getMaxLength_old(listTranscrit):
	length = 0
	for i in listTranscrit:
		if int(i[2]) - int(i[1]) > length:
			length = int(i[2]) - int(i[1])
		else: 
			continue
	return length

'''
Return the maximum length of an a given exon
'''	
def getMaxLength(listTranscrit, dictSequence, lastEnd):
	length = 0
	flag = False
	for i in listTranscrit:		
		#length = int(i[2]) - int(i[1])
		seq = dictSequence[i[0]]
		#print(lastEnd)
		#if i[0] == "ENSTBET00000010850":
			#print(seq)
			#print(lastEnd)
		seq = seq[lastEnd:]

		j = 0
		numberOfNucleotide = 0
		#print(seq, int(i[2]))
		for k in range(len(seq)):
			if numberOfNucleotide == (int(i[2])- int(i[1])):
				break

			nucleotide = seq[k]	

			j +=1 
			if nucleotide == "-":
				pass
			else:
				numberOfNucleotide += 1

		if j > length:			
			length = j
			#print(i[0], length)
	#print("=======================================================", length, lastEnd+length)
	return length, lastEnd+length

'''
Get the sequence of an exon foreach transcript
ID of the transcrit, start and end position of the exon in the transcript, and the length are given
'''	
def getSequence(IDtranscrit, start, end, length, dictSequence, fromSeq, toSeq):
	counterRef = 0 
	finalSequenceRef = ""
	
	if int(end) == 0 and int(start) == 0: 
		
		#print IDtranscrit, finalSequenceRef
		return IDtranscrit, start, end, finalSequenceRef
	else:
		seq = dictSequence.get(IDtranscrit)
		finalSequenceRef = seq[fromSeq:toSeq]
		endRef = int(start) + int(length)
		return IDtranscrit, start, end, finalSequenceRef

		"""
		sequenceRefExon = dictSequence.get(IDtranscrit)
		endRef = int(start) + int(length)
		for q in sequenceRefExon:
			if counterRef < int(endRef): 
				if q == '-' and counterRef < int(start) :
					pass 
				elif q != '-' and counterRef >= int(start): 
					counterRef += 1
					finalSequenceRef += q
				elif q == '-' and counterRef >= int(start):
					counterRef += 1
					finalSequenceRef += q
				else :
					counterRef += 1
		#print IDtranscrit, finalSequenceRef		
		return IDtranscrit, start, endRef, finalSequenceRef
		"""


'''
For a given transcript, return the length and the phase of the exon
'''	
def getPhaseLength(seqTranscrit, start):
	newSeqTranscrit = seqTranscrit.replace("-", "") 
	phase = abs(int(start))%3
	length = len(newSeqTranscrit)
	return phase, length

'''
Foreach exon we have a list with the IDtranscrit, start and end position in the transcript of the exon and the microalignment sequence of the exon for the given transcript
'''	
def sequenceDictionary(macro, dictSequence, sequenceFinal):
	lastEnd = 0
	for n in range(len(macro)):
		templist = []
		#print(macro[n])
		fromSeq = lastEnd
		length, lastEnd = getMaxLength(macro[n], dictSequence, lastEnd)
		toSeq = lastEnd
		#print(length, lastEnd)

		for m in range(len(macro[0])): 
			#print(macro[n][m][0], macro[n][m][1], macro[n][m][2])

			IDtranscrit, start, end, sequenceExon = getSequence(macro[n][m][0], macro[n][m][1], macro[n][m][2], length, dictSequence, fromSeq, toSeq)
			#print(IDtranscrit, start, end)
			#ENSTBET00000010850
			#if IDtranscrit == "ENSAHAT00000002902":
				#print( macro[n][m])				
				#print (IDtranscrit, start, end ,sequenceExon)
				#print(macro[n][m][0], macro[n][m][1], macro[n][m][2], length, fromSeq, toSeq)
				#print("_______________________________________________________________________\t")
			templist.append([IDtranscrit, start, end ,sequenceExon])
		#print(templist)		
		#if length == 1181:
		#	print (macro[n])#, dictSequence, lastEnd)
		#	exit()
		sequenceFinal.append(templist)
	#exit()
'''
Skipping exon event 
'''	
def defineEventExonSkipping(transcrit, start, end):
	return True if int(start) == 0 and int(end) == 0 else False

'''
Only gap event 
'''	
def onlyGapEvent(transcrit, start, end):
	nbGap = 0 
	if int(start) == 0 and int(end) == 0:
		return False
	else:
		for gap in transcrit: 
			if gap == '-':
				nbGap += 1 
	return True if nbGap == len(transcrit) else False
	

'''
3 prime event 
'''	
def defineEvent3(transcrit):
	tprime = 0
	for o in reversed(transcrit):
		if o == '-' :
			tprime += 1
		else:
			break
	return True if tprime > 0 else False

'''
5 prime event 
'''	
def defineEvent5(transcrit):
	cprime = 0
	for o in transcrit:
		if o == '-' :
			cprime += 1
		else:
			break
	return True if cprime > 0 else False

'''
Internal gap <> 30 gap event 
'''	
def defineEventInternalGap(transcrit):
	listGap = []
	positionStart = 0 
	positionEnd = 0
	length = 0
	flagEndGap = True
	
	for nucleotide in transcrit:
		if nucleotide == '-':
			length += 1
			positionEnd = positionStart + 1
		else:
			positionStart += 1
			continue
		positionStart += 1
		
	if int(length) == 0 or int(positionEnd) == 0 or int(positionEnd - length) == 0 or int(positionEnd) == int(len(transcrit)):
		return 0
	else :		
		listGap.append([positionEnd - length, positionEnd, length])
		return listGap
	
'''
No event the exon is preserved
'''	
def defineEventPreserved(transcrit):	
	return False if defineEvent3(transcrit) and defineEvent5(transcrit) or defineEventInternalGap(transcrit) != 0 else True



def main_labelExons(save_path, macroalignment, microalignment):
	dictSequence = {}
	macro = []
	sequenceFinal = []	
	namefile = macroalignment.split("/")[-1]
	namefile = namefile.split("_")[0]
	namefile = namefile + "_labeledExons.out"
	completeName = os.path.join(save_path, namefile)         
	outputfile = open(completeName, "w")	
		
	microalignmentIntoDictionnary(microalignment, dictSequence)
	macroalignmentIntoList(macroalignment, macro)
	sequenceDictionary(macro, dictSequence, sequenceFinal)
	

	numExon = 0
	for exon in sequenceFinal:
		outputfile.write(">" + str(numExon) + "\n")
		for transcrit in exon:
			event = False
			res = []
			outputfile.write(transcrit[0] + ": ")
			phase, length = getPhaseLength(transcrit[3], transcrit[1])
			#if transcrit[0] == "ENSTBET00000010850":
				#print(transcrit[3], transcrit[2], transcrit[1])
			if defineEventExonSkipping(transcrit[3], transcrit[2], transcrit[1]):
				event = True
				res.append("exon skipping")
			if onlyGapEvent(transcrit[3], transcrit[2], transcrit[1]):
				#print((transcrit[3], transcrit[2], transcrit[1]))
				#dans le fichier macro alignment, les positions de debut et fin des cds ne commencent pas a 0 et
				#et ne sont pas consecutives.
				#print("---------------------------------------------------------------------------")
				#exit("ONLY GAP")
				event = True
				res.append("only gap")
			else :
				if defineEvent3(transcrit[3]):
					event = True
					res.append("3event")
				if defineEvent5(transcrit[3]):
					event = True
					res.append("5event")
				if defineEventInternalGap(transcrit[3]) != 0:
					event = True
					l = defineEventInternalGap(transcrit[3])
					if l[0][2] >= 30 :
						res.append(["internal gap >= 30", l[0][0], l[0][2]])
					else: 
						res.append(["internal gap < 30", l[0][0], l[0][2]])

			if not event: 
				res.append("no event")
				
			res.append(length)
			res.append(phase)
			outputfile.write(str(res))
			outputfile.write("\n")
		numExon += 1

	outputfile.close()
	return completeName


def main_labelExons_aux(macroalignment, microalignment):
	path =  sys.path[0]      
	save_path  = path + "/labelExons/"	
	return main_labelExons(save_path, macroalignment, microalignment)

if __name__ == "__main__":    
	path =  sys.path[0]
	parser = build_arg_parser(path)
	arg = parser.parse_args()        
	save_path  = arg.save_path 
	macroalignment = arg.macroalignment
	microalignment = arg.microalignment

	main_labelExons(save_path, macroalignment, microalignment)
