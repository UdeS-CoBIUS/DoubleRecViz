# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

'''

import sys
from ete3 import Tree
from functions import *

def prune_species_tree(file_input_tree,file_taxa_list):
	# input : name_of_input.tree, list_of_taxa_to_keep.txt
	# output : newick string of species tree
	input_tree = get_input_tree(file_input_tree)
	taxa_list = get_taxa_list(file_taxa_list)
	t = Tree(input_tree,format=1)
	t.prune(taxa_list)
	node_names = []
	i = 0
	for node in t.traverse("postorder"):
		if(node.name == ""):
			node.name = "S"+str(i)
			i+=1
		while(node.name in node_names):
			node.name += "1"
		node_names.append(node.name)
			
	s = t.write(format=8)
	s = add_root_nw(s,"S")
	return s
        
def get_taxa_list(file_taxa_list) : 
	taxa_list = []
	try :
		f = open(file_taxa_list, 'r')
		for line in f:
			taxa_list.append(line.split('\n')[0].split(" ")[-1])
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	return taxa_list

