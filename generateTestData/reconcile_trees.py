# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

'''

import os
from functions import *

def reconcile(tree1,tree2, file_mapping):
	t1 = add_content_tree(tree1,file_mapping)
	labeled_tree1 = "labeled.nw"
	labeled_tree2 = "result.nw"
	write_tree(t1,labeled_tree1,'w')
	os.system("treebest sdi -s "+tree2 + " " + labeled_tree1 + " > " + labeled_tree2)
	os.system("rm "+labeled_tree1)
	s = read_nw(labeled_tree2)
	os.system("rm "+labeled_tree2)
	return s
