# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

'''

from Bio import SeqIO
from Bio.Seq import translate
import os
from ete3 import Tree


def add_content_tree(file_input_tree,file_mapping):
	# input : name_of_input.tree, mapping.txt
	# output : newick string of gene tree
	input_tree = get_input_tree(file_input_tree)
	mapping = get_mapping(file_mapping)
	t = Tree(input_tree, format = 1)
	i = 0
	for node in t.traverse("postorder"):
		if(node.name in mapping.keys()):
			node.name = node.name + "_"+mapping[node.name]
	s = t.write(format=8)
	return s 

def add_root_nhx(input_tree,prefix):
	s1,s2 = input_tree.split(")[")
	s = s1 + ")"+prefix+"ROOT["+s2
	return s 

def add_root_nw(input_tree,prefix):
	s1,s2 = input_tree.split(";")
	s = s1 +prefix+"ROOT;"
	return s 

def remove_content_tree_nhx(file_input_tree):
	input_tree = get_input_tree(file_input_tree)
	root_feature = "["+input_tree.split(")[")[1]
	t = Tree(input_tree, format = 1)
	i = 0
	for node in t.traverse("postorder"):
                node.name = node.name.split("_")[0]
	s = t.write(format=8, features=[])
	s = s.split(";")[0]+root_feature
	return s 

def remove_content_tree_nw(file_input_tree):
	input_tree = get_input_tree(file_input_tree)
	t = Tree(input_tree, format = 1)
	i = 0
	for node in t.traverse("postorder"):
                node.name = node.name.split("_")[0]
	s = t.write(format=8)
	return s 

def get_mapping(file_mapping) : 
	mapping = {}
	try :
		f = open(file_mapping, 'r')
		for line in f:
			id1,id2 = line.split('\n')[0].split(" ")
			mapping[id1] = id2
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))
	return mapping

def read_nw(file_input_tree):
	# input : name_of_input.tree
	# output : newick string of gene tree
	s = ""
	try :
		f = open(file_input_tree, 'r')
		for line in f:
			ss = line.split('\n')[0]
			if(len(ss) > 0):
				s += ss
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))
	return s 

def write_tree(string_input_tree,file_input_tree, mode):
	f = open(file_input_tree, mode)
	f.write(string_input_tree+"\n")
	f.close()

def read_recGeneTreeXML(file_input_tree):
	s = ""
	try :
		f = open(file_input_tree, 'r')
		for line in f:
			if("recPhylo" not in line and
                           "confidence" not in line and
                           "branch_length" not in line):
				s += line
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))
	return s 

def read_recTransTreeXML(file_input_tree):
	s = ""
	try :
		f = open(file_input_tree, 'r')
		for line in f:
			if("recPhylo" not in line and
                           "confidence" not in line and
                           "branch_length" not in line and
                           "recGeneTree" not in line):
				if("duplication" in line):
					line = line.replace("duplication", "creation")
				if("speciesLocation" in line):
					line = line.replace("speciesLocation", "genesLocation")
				s += line
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))
	s = "  <recTransTree>\n" + s + "  </recTransTree>\n"
	return s 

def compute_alignment(sequences,file_mapping):
	AAsequences = sequences.split('.fasta')[0] + "AA.fasta"
	AAfile = open(AAsequences,"w")
	for record in SeqIO.parse(sequences, "fasta"):
		id = record.id
		seq = str(translate(str(record.seq)))
		AAfile.write(">"+id+"\n")
		AAfile.write(seq+"\n")        
	AAfile.close()

	AAsequences_aln = AAsequences.split('.fasta')[0] + ".aln"
        
	command = "muscle -in "  + AAsequences + " -out " +  AAsequences_aln + " 2>/dev/null >/dev/null"
	os.system(command)	
	os.system("rm "+ AAsequences)	

	nt_dict = {}
	for record in SeqIO.parse(sequences, "fasta"):
		nt_dict[record.id] = str(record.seq)

	aa_aln_dict = {}
	for record in SeqIO.parse(AAsequences_aln, "fasta"):
		aa_aln_dict[record.id] = str(record.seq)		

	os.system("rm "+ AAsequences_aln)	

	sequences_aln = sequences.split('.fasta')[0] + ".aln"
	mapping = {}
	try :
		f = open(file_mapping, 'r')
		for line in f:
			id1,id2 = line.split('\n')[0].split(" ")
			mapping[id1] = id2
	except IOError as e:
                print("I/O error({0}): {1}".format(e.errno, e.strerror))

	file = open(sequences_aln, "w")
	for id1 in aa_aln_dict.keys():
		back_translate_seq = ""
		aa_aln_seq = aa_aln_dict[id1]
		nt_seq = nt_dict[id1]
		nb_gap = 0
		for pos in range(len(aa_aln_seq)):
			if aa_aln_seq[pos] == "-":
				back_translate_seq += "---"
				nb_gap += 1
			else:
				pos_to_use = pos - nb_gap
				back_translate_seq += nt_seq[3*pos_to_use:3*pos_to_use+3]
		file.write(">"+id1+"_"+mapping[id1]+"\n")
		file.write(back_translate_seq + "\n")

	file.close()
	return sequences_aln

def  get_input_tree(file_input_tree):
	input_tree = ""
	try :
		f = open(file_input_tree, 'r')
		line = f.readline()
		input_tree = line.split('\n')[0]
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	return input_tree
