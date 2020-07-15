# -*- coding: utf-8 -*-
from shutil import copyfile
import argparse

def build_arg_parser():
    parser = argparse.ArgumentParser(description="move gene family files")
    parser.add_argument('-g', '--genefamilyid')
    return parser


def moveFile(gene_familly_id):
	copyfile("SpliceFamAlignMulti/examples/output/" + gene_familly_id + "_initialsource.fasta", "SHT/initialSource/" + gene_familly_id + "_initialsource.fasta")
	copyfile("SpliceFamAlignMulti/examples/output/" + gene_familly_id + "_initialsource2target.txt", "SHT/initialSource/" + gene_familly_id + "_initialsource2target.txt")	
	copyfile("SpliceFamAlignMulti/examples/output/" + gene_familly_id + "_macroalignment.txt", "SHT/macroalignment/" + gene_familly_id + "_macroalignment.txt")
	copyfile("SpliceFamAlignMulti/examples/output/" + gene_familly_id + "_microalignment.fasta", "SHT/microalignment/" + gene_familly_id + "_microalignment.fasta")
	copyfile("SpliceFamAlignMulti/examples/output/" + gene_familly_id + "_initialsource.fasta", "SHT/initialSource/" + gene_familly_id + "_initialsource.fasta")
	

def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	genefamilyid = arg.genefamilyid

	moveFile(genefamilyid)

if __name__ == "__main__":
    main()

