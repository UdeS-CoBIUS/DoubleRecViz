#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import requests
import os
from os import getcwd
import argparse
from ete3 import Tree
import glob
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import random, string
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
#from StringIO import StringIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from copy import deepcopy, copy
import traceback
import matplotlib.pyplot as plt
import pandas as pd

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


def detailGeneTranscript(mappingGeneToTranscriptInit):
	mappingGeneToTranscript = {}
	for k, v in mappingGeneToTranscriptInit.items():
		for seqId in v:
			mappingGeneToTranscript[seqId] = k
	return mappingGeneToTranscript
	
i = 0
for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/clusters_100_0/fuzzyCMeans/*_cluster_1_root.nhx"): 
    i += 1 
    x = str(x)
    x= x.split("/")
    x = x[-1]		
    x = x.split("_")[0]	
    mappingGeneToTranscript = {}
    intial_source = "initialSource/" + x + "_initialsource.fasta"
    intialsource2target = open("initialSource/" + x + "_initialsource2target.txt", "r")
    intial_source_AA = open("initialSource/" + x + "_initialsource_isosel_AA.fasta", "w")
    intial2source_isosel = open("initialSource/" + x + "_initialsource2target_isosel.txt", "w")
    
    for record in SeqIO.parse(intial_source, "fasta"):
        id = record.id
        seq = str(translate(record.seq))
        intial_source_AA.write(">" + id + "\n")
        intial_source_AA.write( seq + "\n\n")        
    intial_source_AA.close()
    
    lines = intialsource2target.readlines()
    for l in lines:
        l.replace("\n", "")
        p = l.split(" ")
        p[1] = p[1].replace("\n", "")
        intial2source_isosel.write(">" + p[0] + "\t" + p[1] + "\n")
        if p[1] in mappingGeneToTranscript.keys():
            mappingGeneToTranscript[p[1]].append(p[0])
        else:
            mappingGeneToTranscript[p[1]] = [p[0]]
            
    intial2source_isosel.close()
    cmd = "IsoSel_source_code2/bin/run_isosel initialSource/" + x + "_initialsource_isosel_AA.fasta -a muscle -b 1 -f initialSource/" + x + "_initialsource2target_isosel.txt --outdir isoselSeq -n " + x
    print(cmd)   
    os.system(cmd)
    
    ids = []
    if os.path.exists("isoselSeq/" + x + "_filtered.fasta"):
        for record in SeqIO.parse("isoselSeq/" + x + "_filtered.fasta", "fasta"):
            id = record.id
            ids.append(id)
            

	    #try:		
	
        mappingGeneToTranscriptDetail = detailGeneTranscript(mappingGeneToTranscript)
        geneToSpecieFIle = open("geneToSpecie/cleanTree/" + x + "_cleanoutput.out", "r")

        mapping_file = open("tmp/" + x + "_mapping_species.txt", "w")
        mapping_transcript_file = open("tmp/" + x + "_mapping_transcript.txt", "w")
        microalignmentFile = "microalignment/" + x + "_microalignment.fasta"

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

	
        seqOfThisCLuster = []
        AASeqOfThisCLuster = []
        for seqId in mappingGeneToTranscriptDetail.keys():
	        sequenceID = seqId + "" + mappingGeneToTranscriptDetail[seqId] + "_" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]]			
	        mapping_transcript_file.write(sequenceID + "\t" + seqId + "\t" + mappingGeneToTranscriptDetail[seqId] +  "\t" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]] +"\n")
        mapping_file.close()
        mapping_transcript_file.close()
        #le programme enregistre uniquement les id des seq qui sont retenus par le master cluster
		
        for seqId in ids:
	        sequence = str(sequences[seqId])
	        sequence = Seq(sequence.replace("-", ""))
	
	        sequenceID = seqId + "" + mappingGeneToTranscriptDetail[seqId] + "_" + speciesOfGene[mappingGeneToTranscriptDetail[seqId]]
	


	        record = SeqRecord(sequence , sequenceID, '', '')
	        seqOfThisCLuster.append(record)

	        AARecord = SeqRecord(translate(sequence) , sequenceID, '', '')
	        AASeqOfThisCLuster.append(AARecord)				
	
	        filename_tmp = "isoselSeq/" + x + ".fasta"		
	        SeqIO.write(seqOfThisCLuster, filename_tmp, "fasta")

	        filename_tmp2 = "isoselSeq/" + x  + "_seqAA.fasta"		
	        SeqIO.write(AASeqOfThisCLuster, filename_tmp2, "fasta")

	        nt_alignment = "isoselSeq/" + x  + "_NT.fasta"
	        aa_alignment = "isoselSeq/" + x  + "_AA.fasta"

	        compute_alignment_muscle(filename_tmp, filename_tmp2, nt_alignment, aa_alignment)



	        sequences_nt = []
	        for record in SeqIO.parse(nt_alignment, "fasta"):
		        sequence = str(record.seq)
		        sequence = Seq(sequence.replace("!", "A"))			
		        sequences_nt.append(SeqRecord(sequence, record.id, '', ''))  

	        SeqIO.write(sequences_nt, nt_alignment, "fasta")


	        command2 = "./treebest best " + nt_alignment + " > " +  "isoselSeq/" + x  + "_root.nhx -f ressources/ensemblSpecies.tree"
	        print (command2)
	        os.system(command2)


	        i +=1

        if i == 1:
            exit()

	    #except Exception, e:
	    #    print(e)
	            
    else:
        pass
























