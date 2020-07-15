# -*- coding: utf-8 -*-
from shutil import copyfile
import re
import argparse
from ete3 import Tree
import glob

def build_arg_parser():
    parser = argparse.ArgumentParser(description="move gene family files")
    parser.add_argument('-g', '--genefamilyid')
    return parser


def extract_trees(gene_family_id):
	file = open("SHT/clusters_50_50/fuzzyCMeans/" + gene_family_id + "_cluster_1_root.nhx", "r")
	tree_str = file.read()
	tree_str = tree_str.replace("\n", "")
	tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
	tree = Tree(tree_str)
	mapping_genes = {}
	mapping_genes_file = open("SHT/tmp/" + gene_family_id + "_mapping_transcript.txt", "r")
	lines = mapping_genes_file.readlines()	
	for line in lines:
		line = line.replace("\n", "")
		line = line.split("\t")
		mapping_genes[line[0]] = [line[1], line[2], line[3]]

	for node in tree.traverse("postorder"):
		if node.is_leaf():
			name = node.name
			node.name =  mapping_genes[name][1]

	tree.write(outfile = "SuperProteinTree/datas/" + gene_family_id + "_genetree.nw", format=9)

	for x in glob.glob("SHT/clusters_50_50/fuzzyCMeans/" + gene_family_id + "_transcrit_*_root.nhx"):			
		copyfile(x, "SuperProteinTree/datas/" + x.split("/")[-1])
	

def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	genefamilyid = arg.genefamilyid
	extract_trees(genefamilyid)

if __name__ == "__main__":
    main()
	

