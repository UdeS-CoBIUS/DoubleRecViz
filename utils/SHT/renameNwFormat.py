#!/usr/bin/python

'''
Author : Marie Degen
Example : python phyml_tree.py /home/local/USHERBROOKE/degm2303/Documents/SpliceGraph/clusters/sequencesClusters


'''

import os
import sys

				
for root, dirs, files in os.walk(sys.argv[1]):  
	for d in dirs:
		for files in os.walk(sys.argv[1] + "/" +d): 
			for i in files:
				for filename in i:
					print filename
					if filename.endswith(".phylip_phyml_tree"):
						newname = filename.split(".")[0] + "nw"
						print newname
						os.rename(filename, newname)
