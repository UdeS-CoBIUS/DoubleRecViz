#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob

def geneToSpecies(clean_original_gene_tree):
	file = open(clean_original_gene_tree, "r")
	geneToSpeciesDict = {}
	lines = file.readlines()
	for line in lines:
		line = line.replace("\n", "")
		line = line.split(" ")
		geneToSpeciesDict[line[1]] = line[0]
	return geneToSpeciesDict



listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
#listes = [[0.5,0.5]]		


for e in  listes:
	weightStructure = e[0]
	weightSequence = e[1]
	sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))

	astral_input_tree_orthogroup = open("Astral/test_data/" + sub_path + "_inputAstralOrthoGroup.tre", "w")
	astral_input_tree_Ensembl = open("Astral/test_data/" + sub_path + "_inputAstralEnsembl.tre", "w")
	astral_input_tree_ReconstructEnsembl = open("Astral/test_data/" + sub_path + "_inputReconstructEnsembl.tre", "w")
	astral_input_tree_IsoSel = open("Astral/test_data/" + sub_path + "_inputIsoSel.tre", "w")


	for root, dirs, files in os.walk("ensemblGeneTree/trees_" + sub_path): 		
			for file in files:					
				trees = open("ensemblGeneTree/trees_" + sub_path + "/" + file, "r")
				lines = trees.readlines()
				listestmp = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
				#listes = [[0.5,0.5]]		

				skip = False				 
				for e in  listestmp:
					weightStructure = e[0]
					weightSequence = e[1]
					sub_pathtmp = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))				
					if os.path.exists("ensemblGeneTree/trees_" + sub_pathtmp + "/" + file) == False:
						skip = True

				if len(lines) != 5 or skip:
					pass
				else:			
					print(file, sub_path)		
					treeFromtreeBest = PhyloTree(lines[0])
					treeFromEnsembl = PhyloTree(lines[2])
					treeFromReconstructEnsembl = PhyloTree(lines[3])
					treeFromIsosel = PhyloTree(lines[4])		
					nb = len(treeFromtreeBest.get_leaf_names())
					if(len(treeFromEnsembl.get_leaf_names()) == nb  and len(treeFromReconstructEnsembl.get_leaf_names()) == nb and len(treeFromIsosel.get_leaf_names()) == nb):

						clean_original_gene_tree = "ensemblGeneTree/cleanTree/" + file.split(".")[0] + "_cleanoutput.out"

						geneToSpeciesDict = geneToSpecies(clean_original_gene_tree)					
						#print(geneToSpeciesDict)
						for node in treeFromtreeBest.traverse("preorder"):
							if node.is_leaf():
								node.name = geneToSpeciesDict[node.name]

						for node in treeFromEnsembl.traverse("preorder"):
							if node.is_leaf():
								node.name = geneToSpeciesDict[node.name]							

						for node in treeFromReconstructEnsembl.traverse("preorder"):
							if node.is_leaf():
								node.name = geneToSpeciesDict[node.name]						

						for node in treeFromIsosel.traverse("preorder"):
							if node.is_leaf():
								node.name = geneToSpeciesDict[node.name]								

						if len(list(set(treeFromReconstructEnsembl.get_leaf_names()))) == len(treeFromReconstructEnsembl.get_leaf_names()):
							astral_input_tree_orthogroup.write(treeFromtreeBest.write(format=9) + "\n")
							astral_input_tree_Ensembl.write(treeFromEnsembl.write(format=9) + "\n")
							astral_input_tree_ReconstructEnsembl.write(treeFromReconstructEnsembl.write(format=9) + "\n")
							astral_input_tree_IsoSel.write(treeFromIsosel.write(format=9) + "\n")

	astral_input_tree_orthogroup.close()
	astral_input_tree_Ensembl.close()
	astral_input_tree_ReconstructEnsembl.close()
	astral_input_tree_IsoSel.close()
	print("java -jar Astral/astral.5.7.3.jar -i Astral/" + "Astral/test_data/" + sub_path + "_inputAstralOrthoGroup.tre " + " -q -o Astral/test_data/" + sub_path + "_SpeciesTreeOrthoGroup.tre -t 16")

	os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + sub_path + "_inputAstralOrthoGroup.tre " + " -o Astral/test_data/" + sub_path + "_SpeciesTreeOrthoGroup.tre")
	os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + sub_path + "_inputAstralEnsembl.tre " + " -o Astral/test_data/" + sub_path + "_SpeciesTreeEnsembl.tre")
	os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + sub_path + "_inputReconstructEnsembl.tre " + " -o Astral/test_data/" + sub_path + "_SpeciesTreeReconstructEnsembl.tre")
	os.system("java -jar Astral/astral.5.7.3.jar -i " + "Astral/test_data/" + sub_path + "_inputIsoSel.tre " + " -o Astral/test_data/" + sub_path + "_SpeciesTreeIsoSel.tre")

