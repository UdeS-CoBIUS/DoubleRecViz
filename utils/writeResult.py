# -*- coding: utf-8 -*-
from shutil import copyfile
import argparse
import os	

def build_arg_parser():
    parser = argparse.ArgumentParser(description="move gene family files")
    parser.add_argument('-id', '--genefamilyid')
    return parser


def writeResults(gene_familly_id):
	try: 
		if not(os.path.isdir):
			os.mkdir("Output/" + gene_familly_id) 
	except OSError as error: 
		print(error) 
	
	copyfile("SuperProteinTree/output/" + gene_familly_id + "_superProteinTree.nw", "Output/" + gene_familly_id + "/" + gene_familly_id + "_superProteinTree.nw")

	copyfile("SuperProteinTree/datas/" + gene_familly_id + "_genetree.nw", "Output/" + gene_familly_id + "/" + gene_familly_id + "_genetree.nw")	
	
	copyfile("SuperProteinTree/datas/" + gene_familly_id + "_speciesTree.nw", "Output/" + gene_familly_id + "/" + gene_familly_id + "_speciesTree.nw")	
	print("\n\n===============================RESULT FILES==============================\n")
	print("Output/" + gene_familly_id + "/" + gene_familly_id + "_superProteinTree.nw\n")
	print("Output/" + gene_familly_id + "/" + gene_familly_id + "_genetree.nw\n")
	print("Output/" + gene_familly_id + "/" + gene_familly_id + "_speciesTree.nw\n")

def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	genefamilyid = arg.genefamilyid

	writeResults(genefamilyid)

if __name__ == "__main__":
    main()
