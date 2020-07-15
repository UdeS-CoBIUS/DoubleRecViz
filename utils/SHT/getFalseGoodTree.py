 #!/usr/bin/python

'''
This script return the in output files the name of the file with good and with wrong newick format. 
This script must be put into the directory where the files to test are.

Arguments
--treefile=[path to file containing one gene tree]

Example usage:
python getFalseGoodTree.py

'''

from ete3 import Tree
import requests, sys
import os


goodtreefile = open("goodtreefile.out", "w+")
falsetreefile = open("falsetreefile.out", "w+")
goodtreeformat = 0 
falsetreeformat = 0

for root, dirs, files in os.walk("."):  
    for filename in files:
	    if filename.endswith(".txt"): 
		fil = open(filename, "r")
		datas = fil.readlines()[0]
		try:
			tree = Tree(str(datas))
			goodtreefile.write(filename)
			goodtreefile.write("\n")
			goodtreeformat += 1
		except:
			falsetreefile.write(filename)
			falsetreefile.write("\n")
			falsetreeformat += 1

falsetreefile.close()
goodtreefile.close()

print str(goodtreeformat) + " has the right newick format and " + str(falsetreeformat) + " has an error in their newick format"	

