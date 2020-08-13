# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

input parameter : transcript_sequence_file.fasta transcript2gene_file.txt name_of_output.tree

'''

import sys

sys.path.insert(1, '../Utils/')

import os
from Bio import SeqIO
from Bio.Seq import translate
from ete3 import Tree
from utils import *

def compute_gene_trees(file_transcript,file_transcript2gene,file_gene2species,file_speciestree):
	# input : transcript_sequence_file.fasta transcript2gene_file.tx
	# output : newick string of two gene trees computed using treebest and phyml
	selected_transcript = select_longest_transcript(file_transcript,file_transcript2gene)
	selected_transcript_aln = compute_alignment(selected_transcript,file_gene2species)
	file_output_tree = "tmp.nw"
        
	treebest_command = "treebest best -f "+file_speciestree+" -o " + file_output_tree + " " + selected_transcript_aln + " 2>/dev/null >/dev/null"
	os.system(treebest_command)
	s_treebest = read_tree_gene(file_output_tree)        
	os.system("rm "+file_output_tree)

	phyml_command = "treebest phyml  -o " + file_output_tree + " " + selected_transcript_aln + " 2>/dev/null >/dev/null"
	os.system(phyml_command)
	file_output_tree_rooted = file_output_tree + "_rooted"
	os.system("treebest root "+ file_output_tree +" > " + file_output_tree_rooted)

	s_phyml = read_tree_gene(file_output_tree_rooted)        
	os.system("rm "+file_output_tree)
	os.system("rm "+file_output_tree_rooted)

	os.system("rm "+selected_transcript)
	os.system("rm "+selected_transcript_aln)
	return s_treebest, s_phyml
        
def  select_longest_transcript(file_transcript,file_transcript2gene):
	alltranscripts = {}
	for record in SeqIO.parse(file_transcript, "fasta"):
		alltranscripts[record.id] = str(record.seq)
                
	gene2transcript = {}
	try :
		f = open(file_transcript2gene, 'r')
		for line in f:
			transcript,gene = line.split('\n')[0].split(" ")
			if(gene in gene2transcript.keys()):
				gene2transcript[gene].append(transcript)
			else:
				gene2transcript[gene] = [transcript]
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))

	selected_transcript = file_transcript.split('.fasta')[0] + "longest.fasta"
	outfile = open(selected_transcript,"w")

	for gene in gene2transcript.keys():
		tr = gene2transcript[gene][0]
		for i in range(len(gene2transcript[gene])):
			tri = gene2transcript[gene][i]
			if(len(alltranscripts[tri]) > len(alltranscripts[tr])):
                            tr = tri
		outfile.write(">"+gene+"\n")
		outfile.write(alltranscripts[tr]+"\n")
	outfile.close()
	return selected_transcript

def compute_alignment(selected_transcript,file_gene2species):
	AAselected_transcript = selected_transcript.split('.fasta')[0] + "AA.fasta"
	AAfile = open(AAselected_transcript,"w")
	for record in SeqIO.parse(selected_transcript, "fasta"):
		gene = record.id
		seq = str(translate(str(record.seq)))
		AAfile.write(">"+gene+"\n")
		AAfile.write(seq+"\n")        
	AAfile.close()

	AAselected_transcript_aln = AAselected_transcript.split('.fasta')[0] + ".aln"
        
	command = "muscle -in "  + AAselected_transcript + " -out " +  AAselected_transcript_aln + " 2>/dev/null >/dev/null"
	os.system(command)	
	os.system("rm "+ AAselected_transcript)	

	nt_dict = {}
	for record in SeqIO.parse(selected_transcript, "fasta"):
		nt_dict[record.id] = str(record.seq)

	aa_aln_dict = {}
	for record in SeqIO.parse(AAselected_transcript_aln, "fasta"):
		aa_aln_dict[record.id] = str(record.seq)		

	os.system("rm "+ AAselected_transcript_aln)	

	selected_transcript_aln = selected_transcript.split('.fasta')[0] + ".aln"
	gene2species = {}
	try :
		f = open(file_gene2species, 'r')
		for line in f:
			gene,species = line.split('\n')[0].split(" ")
			gene2species[gene] = species
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))

	file = open(selected_transcript_aln, "w")
	for gene in aa_aln_dict.keys():
		back_translate_seq = ""
		aa_aln_seq = aa_aln_dict[gene]
		nt_seq = nt_dict[gene]
		nb_gap = 0
		for pos in range(len(aa_aln_seq)):
			if aa_aln_seq[pos] == "-":
				back_translate_seq += "---"
				nb_gap += 1
			else:
				pos_to_use = pos - nb_gap
				back_translate_seq += nt_seq[3*pos_to_use:3*pos_to_use+3]
		file.write(">"+gene+"_"+gene2species[gene]+"\n")
		file.write(back_translate_seq + "\n")

	file.close()
	return selected_transcript_aln

def read_tree_gene(file_output_tree):
	s = read_nw(file_output_tree)
	t = Tree(s)
	i = 0
	for node in t.traverse("postorder"):
		if(node.name == ""):
			node.name = "G"+str(i)
			i += 1
		else:
			node.name = node.name.split("_")[0]
	s = t.write(format=8)
	s = add_root_nw(s,"G")
	return s
