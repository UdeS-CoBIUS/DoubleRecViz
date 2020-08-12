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


def copy_gene_tree_init(file_input_tree):
	# input : name_of_input.tree
	# output : newick string of gene tree
	input_tree = get_input_tree(file_input_tree)
	t = Tree(input_tree)
	i = 0
	for node in t.traverse("postorder"):
		if(node.name == ""):
			node.name = "G"+str(i)
			i+=1
	s = t.write(format=8)
	s = add_root_nw(s,"G")
	return s 

def  get_input_tree(file_input_tree):
	input_tree = ""
	try :
		f = open(file_input_tree, 'r')
		line = f.readline()
		input_tree = line.split('\n')[0]
	except IOError as e:
		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	return input_tree
        

