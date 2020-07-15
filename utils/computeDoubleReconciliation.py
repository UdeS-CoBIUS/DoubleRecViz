# -*- coding: utf-8 -*-
from shutil import copyfile
import argparse
from ete3 import PhyloTree
import copy

def build_arg_parser():
	parser = argparse.ArgumentParser(description="Compute Double Reconciliation")
	parser.add_argument('-id', '--genefamilyid')
	return parser

def reduceLeaveName(tree, car):
	dictSpe = {}
	dictSpeInv = {}
	compt = 10
	for node in tree.traverse("preorder"):
		if node.is_leaf():
			var = car + str(compt)
			compt = compt + 1
			dictSpe[node.name] = var
			dictSpeInv[var] = node.name
			node.name = var
	return tree, dictSpe, dictSpeInv

def pruneSpeciesTree(gene_familly_id):
	gene_tree_f = open("SuperProteinTree/datas/" + gene_familly_id + "_genetree.nw", "r")
	gene_tree_content = gene_tree_f.readlines()[0]
	gene_tree = PhyloTree(gene_tree_content)
	gene_leaves = gene_tree.get_leaf_names()

	transcrit_tree_f = open("SuperProteinTree/output/" + gene_familly_id + "_superProteinTree.nw", "r")
	transcrit_content = transcrit_tree_f.readlines()[0]
	transcrit_tree = PhyloTree(transcrit_content)
	
	protein2gene = {}
	gene2specie = {}
	specieslist = []
	trans2gene2species = open("SHT/tmp/" + gene_familly_id + "_mapping_transcript.txt", "r")
	lines = trans2gene2species.readlines()
	for l in lines:	
		l = l.replace("\n", "")
		p = l.split("\t")		
		if p[2] in gene_leaves:
			protein2gene[p[1]] = p[2]
			gene2specie[p[2]] = p[3]
			specieslist.append(p[3])
	speciesTree_f = open("SHT/ressources/ensemblSpecies.tree", "r")
	speciesTree_content = speciesTree_f.readlines()[0]
	speciesTree_content = speciesTree_content.replace("*", "")
	speciesTree = PhyloTree(speciesTree_content, quoted_node_names=False, format=1)
	speciesTree.prune(specieslist)
	file = open("SuperProteinTree/datas/" + gene_familly_id + "_speciesTree.nw","w")
	file.write(speciesTree.write(format=9))
	transcrit_tree_label = copy.deepcopy(transcrit_tree)
	gene_tree_label = copy.deepcopy(gene_tree)
	species_tree_label = copy.deepcopy(speciesTree)
	
	transcrit_tree_label, dictTrans, dictTransInv = reduceLeaveName(copy.deepcopy(transcrit_tree_label), "T")
	gene_tree_newname, dictGene, dictGeneInv = reduceLeaveName(copy.deepcopy(gene_tree_label), "G")
	species_tree_new_name, dictSpecies, dictSpeciesInv = reduceLeaveName(species_tree_label, "S")
	
	for node in transcrit_tree_label.traverse("postorder"):
		if node.is_leaf():
			node.name =  dictGene[protein2gene[dictTransInv[node.name]]] + "_" + node.name

	gene_tree_newname_label = copy.deepcopy(gene_tree_newname)
	for node in gene_tree_newname.traverse("postorder"):	
		if node.is_leaf():
			node.name = dictSpecies[gene2specie[dictGeneInv[node.name]]]  + "_" + node.name

	return dictTrans, dictTransInv, dictGene, dictGeneInv, dictSpecies, dictSpeciesInv, transcrit_tree, transcrit_tree_label, gene_tree_newname_label,  gene_tree, gene_tree_newname, speciesTree, species_tree_new_name

def doubleReconciliation(dictTrans, dictTransInv, dictGene, dictGeneInv, dictSpecies, dictSpeciesInv,transcrit_tree, transcrit_tree_label, gene_tree_newname_label, gene_tree, gene_tree_newname, speciesTree, species_tree_new_name, genefamilyid):
	recon_tree_transcript_gene, events_transcript_gene = transcrit_tree_label.reconcile(gene_tree_newname_label)

	for node in recon_tree_transcript_gene:
		if node.is_leaf():
			p = node.name.split("_")
			if len(p) == 1:
				node.name = dictGeneInv[p[0]]
			elif len(p) == 2:
				node.name = dictGeneInv[p[0]] + "_" + dictTransInv[p[1]]
	print(recon_tree_transcript_gene)
	file = open("Output/" + genefamilyid + "/" + genefamilyid + "_reconciliation_transcripts_genes.nw","w")
	file.write(recon_tree_transcript_gene.write(format=9))

	recon_tree_gene_species, events_gene_species = gene_tree_newname.reconcile(species_tree_new_name)
	for node in recon_tree_gene_species:
		if node.is_leaf():
			p = node.name.split("_")
			if len(p) == 1:
				node.name = dictSpeciesInv[p[0]]
			elif len(p) == 2:
				node.name = dictSpeciesInv[p[0]] + "_" + dictGeneInv[p[1]]
	print(recon_tree_gene_species)
	file = open("Output/" + genefamilyid + "/" + genefamilyid + "_reconciliation_genes_species.nw","w")
	file.write(recon_tree_gene_species.write(format=9))


def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	genefamilyid = arg.genefamilyid

	dictTrans, dictTransInv, dictGene, dictGeneInv, dictSpecies, dictSpeciesInv, transcrit_tree, transcrit_tree_label, gene_tree_newname_label,  gene_tree, gene_tree_newname, speciesTree, species_tree_new_name = pruneSpeciesTree(genefamilyid)

	doubleReconciliation(dictTrans, dictTransInv, dictGene, dictGeneInv, dictSpecies, dictSpeciesInv,transcrit_tree, transcrit_tree_label, gene_tree_newname_label, gene_tree, gene_tree_newname, speciesTree, species_tree_new_name, genefamilyid)

if __name__ == "__main__":
	main()

