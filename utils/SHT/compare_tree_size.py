 #!/usr/bin/python
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
from Bio import Seq
import random, string
from Bio.Align.Applications import MafftCommandline
#from StringIO import StringIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from copy import deepcopy, copy
import traceback



values = []
for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/trees_20_80/*.nw"): 		
	x2 =x
	file = open(x, "r")
	x = str(x)
	x= x.split("/")
	x = x[-1]		
	x = x.split(".")[0]

	lines = file.readlines()
	tree_nw = lines[0]
	tree = Tree(tree_nw)
	nb_leaves = len(tree.get_leaves())

	file2 = open("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/parseTree/" + x + "_output.out", "r")
	lines2 = file2.readlines()
	nb_leaves2 = len(lines2)
	
	if nb_leaves/nb_leaves2 <=1:
		values.append(round(nb_leaves/nb_leaves2, 2))

print(values)