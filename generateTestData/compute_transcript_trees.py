# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke


'''

import sys

sys.path.insert(1, '../Utils/')

import os
from Bio import SeqIO
from Bio.Seq import translate
from ete3 import Tree
from utils import *

def compute_transcript_trees(file_transcript,file_transcript2gene,file_genetree_treebest, file_genetree_phyml):
	transcript_aln = compute_alignment(file_transcript,file_transcript2gene)
	file_output_tree = "tmp.nw"
        
	treebest_best_treebest_command = "treebest best -f "+file_genetree_treebest+" -o " + file_output_tree + " " + transcript_aln + " 2>/dev/null >/dev/null"
	os.system(treebest_best_treebest_command)
	s_treebest_best_treebest = read_tree_transcript(file_output_tree)        
	os.system("rm "+file_output_tree)

	treebest_best_phyml_command = "treebest best -f "+file_genetree_phyml+" -o " + file_output_tree + " " + transcript_aln + " 2>/dev/null >/dev/null"
	os.system(treebest_best_phyml_command)
	s_treebest_best_phyml = read_tree_transcript(file_output_tree)        
	os.system("rm "+file_output_tree)

	phyml_command = "treebest phyml  -o " + file_output_tree + " " + transcript_aln + " 2>/dev/null >/dev/null"
	os.system(phyml_command)
	file_output_tree_rooted = file_output_tree + "_rooted"
	os.system("treebest root "+ file_output_tree +" > " + file_output_tree_rooted)

	s_phyml = read_tree_transcript(file_output_tree_rooted)        
	os.system("rm "+file_output_tree)
	os.system("rm "+file_output_tree_rooted)

	os.system("rm "+transcript_aln)
	return s_treebest_best_treebest, s_treebest_best_phyml, s_phyml
        

def read_tree_transcript(file_output_tree):
	s = read_nw(file_output_tree)
	t = Tree(s)
	i = 0
	for node in t.traverse("postorder"):
		if(node.name == ""):
			node.name = "T"+str(i)
			i += 1
		else:
			node.name = node.name.split("_")[0]
	s = t.write(format=8)
	s = add_root_nw(s,"T")
	return s
