#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
import matplotlib.pyplot as plt
import pandas as pd


def geneToSpecies(clean_original_gene_tree):
	file = open(clean_original_gene_tree, "r")
	geneToSpeciesDict = {}
	lines = file.readlines()
	for line in lines:
		line = line.replace("\n", "")
		line = line.split(" ")
		geneToSpeciesDict[line[1]] = line[0]
	return geneToSpeciesDict

def main_ensenmblTree():	
	for root, dirs, files in os.walk("ensemblGeneTree/geneTree/"):  
		for filename in files:			
			print(filename.split(".")[0])				
			if filename.split(".")[0] not in ["ENSGT00940000166644", "ENSGT00960000192966","ENSGT00940000155807","ENSGT00940000181851","ENSGT00390000004040","ENSGT00390000018635","ENSGT00390000003510","ENSGT00960000191882","ENSGT00390000006370","ENSGT00910000147141","ENSGT00390000018064","ENSGT00390000018351","ENSGT00940000166090","ENSGT00950000182820","ENSGT00940000175222","ENSGT00940000163251","ENSGT00940000159174","ENSGT00530000069597","ENSGT00940000181295","ENSGT00940000157564","ENSGT00960000188289","ENSGT00940000162271","ENSGT00910000146738","ENSGT00960000192520","ENSGT00920000149541","ENSGT00940000155400","ENSGT00910000146926","ENSGT00490000043356","ENSGT00940000176084","ENSGT00940000158217","ENSGT00950000186036","ENSGT00940000166868","ENSGT00940000176102","ENSGT00900000141767","ENSGT00960000189062","ENSGT00950000183551","ENSGT00660000097260","ENSGT00950000186111","ENSGT00960000189633","ENSGT00940000179287","ENSGT00940000167924","ENSGT00950000183099","ENSGT00940000159954","ENSGT00500000045767","ENSGT00950000182910","ENSGT00940000179048","ENSGT00940000171707","ENSGT00660000097387","ENSGT00940000155162","ENSGT00940000164774","ENSGT00940000170265","ENSGT00940000182398","ENSGT00940000162793","ENSGT00940000157805","ENSGT00940000176603","ENSGT00960000191223","ENSGT00960000189018","ENSGT00940000159745","ENSGT00900000143637","ENSGT00940000161903","ENSGT00940000160916","ENSGT00390000002871","ENSGT00910000147645","ENSGT00940000167576","ENSGT00940000170751","ENSGT00940000174048","ENSGT00960000190906","ENSGT00940000178523","ENSGT00940000154628","ENSGT00950000185140","ENSGT00940000163885","ENSGT00940000176406","ENSGT00390000006103","ENSGT00940000159651","ENSGT00940000168920","ENSGT00650000094826","ENSGT00940000169667","ENSGT00910000147974","ENSGT00940000172138","ENSGT00940000162105","ENSGT00940000173382","ENSGT00940000171370","ENSGT00940000173221","ENSGT00940000154956","ENSGT00960000191527","ENSGT00950000185085","ENSGT00940000177785","ENSGT00390000008744","ENSGT00600000085077","ENSGT00960000190203","ENSGT00920000150796","ENSGT00940000173843","ENSGT00390000008939","ENSGT00390000014922","ENSGT00940000160433","ENSGT00940000178205","ENSGT00940000179125","ENSGT00940000179005","ENSGT00940000177095","ENSGT00940000160503","ENSGT00940000174400","ENSGT00390000005224","ENSGT00730000111186","ENSGT00730000111283","ENSGT00390000003337","ENSGT00930000152474","ENSGT00510000049503","ENSGT00940000172622","ENSGT00940000158393","ENSGT00940000166003","ENSGT00650000094861","ENSGT00940000157624","ENSGT00910000146845","ENSGT00960000189229","ENSGT00940000179482","ENSGT00940000179227","ENSGT00940000179803","ENSGT00390000016470","ENSGT00940000176747","ENSGT00940000171603","ENSGT00940000164309","ENSGT00390000001810","ENSGT00940000162636","ENSGT00940000162567","ENSGT00530000066353","ENSGT00390000005412","ENSGT00940000162938","ENSGT00390000006247","ENSGT00940000157987","ENSGT00530000063787","ENSGT00940000172519","ENSGT00660000097391","ENSGT00940000164532","ENSGT00910000148212","ENSGT00390000007890","ENSGT00960000188737","ENSGT00950000183793","ENSGT00530000066781","ENSGT00390000007336","ENSGT00940000179611","ENSGT00940000163053","ENSGT00940000173894","ENSGT00940000180784","ENSGT00640000091497","ENSGT00940000166158","ENSGT00940000153289","ENSGT00660000097235","ENSGT00940000158659","ENSGT00530000066831","ENSGT00940000172532","ENSGT00950000185074","ENSGT00940000162518","ENSGT00920000149868","ENSGT00940000162992","ENSGT00940000175713","ENSGT00940000172369","ENSGT00940000166522","ENSGT00960000192424","ENSGT00940000174692","ENSGT00910000147592","ENSGT00940000169009","ENSGT00940000181708","ENSGT00390000010980","ENSGT00940000178623","ENSGT00940000153460","ENSGT00940000181015","ENSGT00390000000279","ENSGT00960000191439","ENSGT00940000174953","ENSGT00910000148066","ENSGT00940000166723","ENSGT00950000185275","ENSGT00390000005420","ENSGT00950000182735","ENSGT00940000156762","ENSGT00510000051622","ENSGT00940000177487","ENSGT00960000191914","ENSGT00860000135934","ENSGT00940000160391","ENSGT00740000117107","ENSGT00940000177470","ENSGT00440000038370","ENSGT00720000109834","ENSGT00940000177577","ENSGT00940000167230","ENSGT00940000177296","ENSGT00390000014910","ENSGT00940000158566","ENSGT00730000113640","ENSGT00960000192614","ENSGT00960000186821","ENSGT00530000065200","ENSGT00940000154344","ENSGT00940000153414","ENSGT00940000165652","ENSGT00940000182558","ENSGT00940000154911","ENSGT00940000155790","ENSGT00940000176648","ENSGT00940000167812","ENSGT00940000181005","ENSGT00940000170036","ENSGT00940000175205","ENSGT00550000075017","ENSGT00950000183248","ENSGT00390000010309","ENSGT00940000166291","ENSGT00950000186178","ENSGT00940000166439","ENSGT00940000156910","ENSGT00940000161331","ENSGT00940000172885","ENSGT00390000010512","ENSGT00950000184715","ENSGT00540000073666","ENSGT00540000073630","ENSGT00390000015804","ENSGT00940000160004","ENSGT00940000173365","ENSGT00940000164045","ENSGT00390000015386","ENSGT00940000179263","ENSGT00910000147763","ENSGT00460000042247","ENSGT00940000165172","ENSGT00940000162132","ENSGT00940000171977","ENSGT00940000164795","ENSGT00390000016811","ENSGT00960000187298","ENSGT00960000190698","ENSGT00390000004589","ENSGT00940000178467","ENSGT00940000175878","ENSGT00940000155347","ENSGT00940000166277","ENSGT00940000156022","ENSGT00940000167593","ENSGT00940000156549","ENSGT00940000171073","ENSGT00940000158156","ENSGT00940000177889","ENSGT00940000160474","ENSGT00940000161652","ENSGT00390000004165","ENSGT00940000157535","ENSGT00930000152750","ENSGT00390000004541","ENSGT00910000148642","ENSGT00940000165675","ENSGT00920000150635","ENSGT00940000159452","ENSGT00910000147183","ENSGT00940000177245","ENSGT00940000168457","ENSGT00940000162164","ENSGT00940000170905","ENSGT00940000170683","ENSGT00920000149145","ENSGT00950000183776","ENSGT00940000178325","ENSGT00940000178526","ENSGT00940000172263","ENSGT00940000158467","ENSGT00940000159081","ENSGT00940000157050","ENSGT00940000158220","ENSGT00940000159129","ENSGT00940000157195","ENSGT00910000148153","ENSGT00910000146522","ENSGT00940000158387","ENSGT00940000155142","ENSGT00940000155294","ENSGT00940000178658","ENSGT00940000164674","ENSGT00390000008536","ENSGT00940000169019","ENSGT00390000002210","ENSGT00940000178224","ENSGT00510000048509","ENSGT00950000185451","ENSGT00940000162542","ENSGT00940000155474","ENSGT00940000169631","ENSGT00940000177500","ENSGT00390000012605","ENSGT00940000155727","ENSGT00910000147055","ENSGT00940000167505","ENSGT00940000153472","ENSGT00940000176997","ENSGT00940000163331","ENSGT00940000161512","ENSGT00940000153692","ENSGT00700000106116","ENSGT00740000115816","ENSGT00550000074935","ENSGT00930000151070","ENSGT00950000183203","ENSGT00860000135889","ENSGT00940000157176","ENSGT00940000163480","ENSGT00940000159800","ENSGT00390000013347","ENSGT00940000181606","ENSGT00940000178764","ENSGT00940000169832","ENSGT00910000144335","ENSGT00950000183312","ENSGT00940000159434","ENSGT00390000003246","ENSGT00390000014792","ENSGT00940000155818","ENSGT00950000185306","ENSGT00940000164366","ENSGT00940000155037","ENSGT00940000176325","ENSGT00390000004717","ENSGT00940000178277","ENSGT00940000157526","ENSGT00390000016298","ENSGT00940000171125","ENSGT00900000143807","ENSGT00510000048991","ENSGT00940000158046","ENSGT00940000169600","ENSGT00940000178375","ENSGT00530000066614","ENSGT00940000169068","ENSGT00940000157707","ENSGT00940000161450","ENSGT00940000156413","ENSGT00940000163940","ENSGT00940000167355","ENSGT00940000158754","ENSGT00940000172014","ENSGT00390000009864","ENSGT00940000162075","ENSGT00940000179368","ENSGT00390000008931","ENSGT00410000029004","ENSGT00390000003571","ENSGT00930000152993","ENSGT00940000159126","ENSGT00880000138575","ENSGT00940000157867","ENSGT00390000006227","ENSGT00900000143393","ENSGT00960000190786","ENSGT00910000147156","ENSGT00510000049083","ENSGT00390000014476","ENSGT00940000178791","ENSGT00940000160970","ENSGT00950000183038","ENSGT00940000159139","ENSGT00950000182863","ENSGT00390000008796","ENSGT00940000168285","ENSGT00940000155823","ENSGT00940000156453","ENSGT00940000165005","ENSGT00960000190855","ENSGT00940000176078","ENSGT00940000162816","ENSGT00940000155217","ENSGT00940000159007","ENSGT00940000157815","ENSGT00390000007427","ENSGT00910000147528","ENSGT00940000156370","ENSGT00940000180465","ENSGT00940000155862","ENSGT00910000146517","ENSGT00940000162522","ENSGT00400000024336","ENSGT00910000148173","ENSGT00940000159695","ENSGT00940000171677","ENSGT00940000159374","ENSGT00940000157459","ENSGT00940000160546","ENSGT00950000186428","ENSGT00910000147567","ENSGT00940000160730","ENSGT00950000183533","ENSGT00940000162734","ENSGT00940000177143","ENSGT00900000143818","ENSGT00940000164615","ENSGT00900000143672","ENSGT00390000013953","ENSGT00390000001511","ENSGT00940000156265","ENSGT00940000155167","ENSGT00940000156791","ENSGT00530000069601","ENSGT00530000067998","ENSGT00940000166658","ENSGT00940000160978","ENSGT00940000162197","ENSGT00940000158783","ENSGT00950000185540","ENSGT00940000176957","ENSGT00940000157071","ENSGT00940000172611","ENSGT00600000085600","ENSGT00660000096058","ENSGT00940000158642","ENSGT00910000146899","ENSGT00940000170385","ENSGT00950000183132","ENSGT00390000013286","ENSGT00940000160763","ENSGT00960000189321","ENSGT00530000063255","ENSGT00940000166644","ENSGT00940000158899","ENSGT00390000003308","ENSGT00940000177899","ENSGT00940000168794","ENSGT00940000161737","ENSGT00940000163329","ENSGT00940000159656","ENSGT00940000169324","ENSGT00940000158362","ENSGT00910000147321","ENSGT00390000003628","ENSGT00960000186606","ENSGT00940000178019","ENSGT00390000014590","ENSGT00860000135092","ENSGT00420000029792","ENSGT00940000158667","ENSGT00940000154730","ENSGT00690000103118","ENSGT00940000162882","ENSGT00940000168313","ENSGT00940000168688","ENSGT00940000179298","ENSGT00720000108834","ENSGT00940000178817","ENSGT00390000010017","ENSGT00390000014331","ENSGT00960000190119","ENSGT00940000163213","ENSGT00940000169756","ENSGT00940000178120","ENSGT00940000158027","ENSGT00940000156489","ENSGT00940000165548","ENSGT00770000121806","ENSGT00940000161156","ENSGT00940000164208","ENSGT00950000183083","ENSGT00940000178799","ENSGT00390000013662","ENSGT00940000164058","ENSGT00940000178042","ENSGT00940000180073","ENSGT00390000002329","ENSGT00940000175300","ENSGT00960000193250","ENSGT00960000192503","ENSGT00960000192228","ENSGT00960000191962","ENSGT00960000191874","ENSGT00960000190584","ENSGT00960000190528","ENSGT00960000189642","ENSGT00960000189089","ENSGT00960000188411","ENSGT00960000187756","ENSGT00950000186274","ENSGT00950000185478","ENSGT00940000182266","ENSGT00940000181424","ENSGT00940000181178","ENSGT00940000180960","ENSGT00940000180934","ENSGT00940000180699","ENSGT00940000179907","ENSGT00940000177773","ENSGT00940000176921","ENSGT00940000174991","ENSGT00940000168127","ENSGT00940000166612","ENSGT00940000163467","ENSGT00940000161544","ENSGT00940000159492","ENSGT00940000157965","ENSGT00940000157919","ENSGT00940000157441","ENSGT00910000148697","ENSGT00900000142991","ENSGT00770000121700","ENSGT00740000116816","ENSGT00730000111779","ENSGT00730000111022","ENSGT00660000097430","ENSGT00650000094935","ENSGT00390000014500","ENSGT00390000008296","ENSGT00390000007041","ENSGT00960000192352","ENSGT00960000190547","ENSGT00960000190521","ENSGT00960000187700","ENSGT00960000187527","ENSGT00960000187439","ENSGT00950000186234","ENSGT00940000182586","ENSGT00940000181310","ENSGT00940000177692","ENSGT00940000177117","ENSGT00940000172915","ENSGT00940000172840","ENSGT00940000171981","ENSGT00940000169608","ENSGT00940000167185","ENSGT00940000164430","ENSGT00940000160714","ENSGT00930000152259","ENSGT00800000125338","ENSGT00730000114600","ENSGT00720000108860","ENSGT00690000103799","ENSGT00680000100706","ENSGT00660000097470","ENSGT00530000068923","ENSGT00390000014175","ENSGT00390000012700","ENSGT00390000007479"]:
				continue	
				
			print(filename.split(".")[0])			
			original_gene_tree = root + filename
			
			clean_original_gene_tree = "ensemblGeneTree/cleanTree/" + filename.split(".")[0] + "_cleanoutput.out"

			geneToSpeciesDict = geneToSpecies(clean_original_gene_tree)
			
			#print(geneToSpeciesDict)
			file_tree = open(original_gene_tree, "r")

			tree_str = file_tree.readlines()[0]

			tree_str = tree_str.replace("/", "")

			tree_str = tree_str.replace("(CL57BL6)", "CL57BL6")

			tree = Tree(tree_str)
			tree2 = Tree(tree_str)
			#print(geneToSpeciesDict.keys())
			for node in tree.traverse("postorder"):
				if node.is_leaf():
					name = node.name.split("__")
					node.name = name[0] + "_" + geneToSpeciesDict[name[0]]

			tree.write(outfile = "ensemblGeneTree/treeForComparaisonEnsembl/" +filename.split(".")[0] + ".nw", format=9)

			for node in tree2.traverse("postorder"):
				if node.is_leaf():
					name = node.name.split("__")
					node.name = name[0]

			tree2.write(outfile = "ensemblGeneTree/treeForComparaisonEnsembl2/" +filename.split(".")[0] + ".nw", format=9)			
			#tree.show()

			#exit()
			#return parseOneFile(save_path, root, filename)



def main_corrected_ensemblTree():	
	error_familly = open("ensemblGeneTree/error_familly.txt", "a")
	cmpt = 0
	listes = [[1,0]]
	#random.shuffle(listes)
	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
		for x in glob.glob("clusters_" + sub_path + "/fuzzyCMeans/*_cluster_1_root.nhx"): 
		#for x in glob.glob("clusters_" + sub_path + "/fuzzyCMeans/ENSGT00390000016811_cluster_1_root.nhx"): 
			try:
				print(x)
				cmpt += 1

		
				mapping_familly = open("ensemblGeneTree/missing_familly.txt", "a")			
				x2 =x

				file = open(x, "r")
				x = str(x)
				x= x.split("/")
				x = x[-1]		
				x = x.split("_")[0]	

				trees_build_by_nj = []
				trees_build_by_orthoGroup = []
				if os.path.exists("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy._seq_OrthoGroup_consel") == False or os.path.exists("ensemblGeneTree/seqFastaUsedByEnsembl/" + x + ".fasta") == False or os.stat("clusters_" + sub_path + "/fuzzyCMeans/" + x + "_nj_" + "1_root.nhx").st_size == 0 or  os.stat("clusters_" + sub_path + "/fuzzyCMeans/" + x + "_cluster_" + "1_root.nhx").st_size == 0:
					pass 
				else:
						
					i = 1
					while os.path.exists("clusters_" + sub_path + "/fuzzyCMeans/" + x + "_nj_" + str(i) + "_root.nhx") and  os.path.exists("clusters_" + sub_path + "/fuzzyCMeans/" + x + "_cluster_" + str(i) + "_root.nhx"):
						trees_build_by_nj.append("clusters_" + sub_path + "/fuzzyCMeans/" + x + "_nj_" + str(i) + "_root.nhx")		
						trees_build_by_orthoGroup.append("clusters_" + sub_path + "/fuzzyCMeans/" + x + "_cluster_" + str(i) + "_root.nhx")
						break
						i += 1
					print("trees_build_by_orthoGroup", trees_build_by_orthoGroup)

					sequencesUsedByOrthoGroup = {}
					for record in SeqIO.parse("initialSource/" + x + "_initialsource.fasta", "fasta"):
						sequencesUsedByOrthoGroup[record.id] = str(record.seq)

					sequencesUsedByEnsemblGeneTree = {}
					for record in SeqIO.parse("ensemblGeneTree/seqFastaUsedByEnsembl/" + x + ".fasta", "fasta"):
						sequencesUsedByEnsemblGeneTree[record.id] = str(record.seq)


					mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
					mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")
					if not (os.path.exists("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt")):
						mapping_familly.write(x + "\n")
						mapping_familly.close()
						continue
					mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")


					mapping_species = {}
					mapping_genes = {}
					geneToseq = {}
					transciptToGene = {}
					mapping_transcript_gene_protein_from_emsenbl = {}
					ensembl_leaves = []
					lines = mapping_transcript_gene_protein_from_emsenbl_file.readlines()
					for line in lines:
						line = line.replace("\n", "")
						parts = line.split("\t")
						if len(parts) == 4:
							mapping_transcript_gene_protein_from_emsenbl[parts[1]] = [parts[0], parts[3], parts[2]]
							ensembl_leaves.append(parts[1])


					lines = mapping_species_file.readlines()	

					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						#print(line)
						mapping_species[line[0]] = line[1]

					lines = mapping_genes_file.readlines()	
					#print(lines)
					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						mapping_genes[line[0]] = [line[1], line[2], line[3]]
						transciptToGene[line[1]] = line[2]


					trees_build_by_nj_nw = []
					trees_build_by_orthoGroup_nw = []

					for tree_st in trees_build_by_nj:
						file = open(tree_st, "r")				
						tree_str = file.read()
						tree_str = tree_str.replace("\n", "")
						tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
						trees_build_by_nj_nw.append(Tree(tree_str))

					for tree_st in trees_build_by_orthoGroup:
						file = open(tree_st, "r")
						tree_str = file.read()
						tree_str = tree_str.replace("\n", "")
						tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
						#print(tree_str)
						trees_build_by_orthoGroup_nw.append(Tree(tree_str))

					#ensembl_tree = Tree("ensemblGeneTree/treeForComparaisonEnsembl2/" + x + ".nw")
					#ensembl_leaves = ensembl_tree.get_leaf_names()

					#Sequence file containing cds choosen among all CDS by orthoGroup approach
					seqFileUsedByOrthoGroup = open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByOrthoGroup.fasta", "w")

					#Sequence file containing the longest CDS of each gene by EnsemblGene Tree
					seqFileUsedByEnsemblGeneTree = open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByEnsemblGeneTree.fasta", "w")
		
					#Here I select just one tree to extract the set of used gene Ids.
					tree = trees_build_by_orthoGroup_nw[0]

					leavesNames = []
					selectedTranscriptUsedByOrthoGroup = []
					for node in tree.traverse("postorder"):
						if node.is_leaf():
							name = node.name
							name =  mapping_genes[name][1]
							if name in ensembl_leaves:
								leavesNames.append(name)
								selectedTranscriptUsedByOrthoGroup.append(mapping_genes[node.name][0])

					#print(ensembl_leaves)
					#print(leavesNames)
					usedGeneIds = []
					for name in transciptToGene.keys():

							geneName = transciptToGene[name]
				
							if name in selectedTranscriptUsedByOrthoGroup:

								seqFileUsedByOrthoGroup.write(">" + geneName + "\n")
								seqFileUsedByOrthoGroup.write(translate(sequencesUsedByOrthoGroup[name]) +  "\n")
								usedGeneIds.append(geneName)

					seqFileUsedByOrthoGroup.close()

					selectedTranscriptUsedByEnsemblGeneTree = []
					if geneName not in mapping_transcript_gene_protein_from_emsenbl.keys():
						continue				
					for geneName in usedGeneIds:
						seqFileUsedByEnsemblGeneTree.write(">" + geneName + "\n")
						proteinIdUsebForThisGeneID =  mapping_transcript_gene_protein_from_emsenbl[geneName][0]
						transcriptIdUsebForThisGeneID =  mapping_transcript_gene_protein_from_emsenbl[geneName][1]
						selectedTranscriptUsedByEnsemblGeneTree.append(transcriptIdUsebForThisGeneID)
						sequences_nt_for_this_proteinID = sequencesUsedByEnsemblGeneTree[proteinIdUsebForThisGeneID]				
						seqFileUsedByEnsemblGeneTree.write(translate(sequences_nt_for_this_proteinID) +  "\n")					


					seqFileUsedByEnsemblGeneTree.close()

					sharedTranscripts = [transcript for transcript in selectedTranscriptUsedByOrthoGroup if transcript in selectedTranscriptUsedByEnsemblGeneTree] 
					notSharedTranscripts1 = [transcript for transcript in selectedTranscriptUsedByOrthoGroup if transcript not in selectedTranscriptUsedByEnsemblGeneTree] 
					notSharedTranscripts2 = [transcript for transcript in selectedTranscriptUsedByEnsemblGeneTree if transcript not in selectedTranscriptUsedByOrthoGroup] 
					totalTranscripts = notSharedTranscripts1 + notSharedTranscripts2 + sharedTranscripts

					stat_file = open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_stats.txt", "w")
					#stat_file.write(str((100.0*len(sharedTranscripts))/len(totalTranscripts)) + "\n")
					#stat_file.write("-".join(sharedTranscripts) +  "\t" + "-".join(notSharedTranscripts1) + "\t" + "-".join(notSharedTranscripts2) + "\t" + "-".join(totalTranscripts))
					stat_file.write(str(len(sharedTranscripts)) +  "\t"+ str(round((100.0*len(sharedTranscripts))/len(selectedTranscriptUsedByOrthoGroup),2)))
					stat_file.close()

					corrected_trees_by_ortho_group = []

					for tree_tmp in trees_build_by_orthoGroup:
						file = open(tree_tmp, "r")
						tree_str = file.read()
						tree_str = tree_str.replace("\n", "")
						tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
						tree = Tree(tree_str)

						parts= tree_tmp.split("/")
						tree_name = parts[-1]		
						tree_name_id = tree_name.split(".nhx")[0]

						#tree.show()
						#gene tree build by Esaie's approach
						tree.write(outfile = "ensemblGeneTree/treeForComparaisonCorrectedEnsembl/" + tree_name_id + ".nw", format=9)

			
						for node in tree.traverse("postorder"):
							if node.is_leaf():				
								name = node.name
								node.name = mapping_genes[name][1]
								if node.name in ensembl_leaves:
									pass
								else:
									node.delete()

						#tree.show()
						tree.write(outfile = "ensemblGeneTree/treeForComparaisonCorrectedEnsembl2/" + tree_name_id + ".nw", format=9)
			
						corrected_trees_by_ortho_group.append(tree.write(format=9))

						leaves_Of_corrected_tree = tree.get_leaf_names()


					corrected_trees_by_ortho_group_nj = []

					for tree_tmp in trees_build_by_nj:
						file = open(tree_tmp, "r")
						tree_str = file.read()
						tree_str = tree_str.replace("\n", "")
						tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
						tree = Tree(tree_str)

						parts= tree_tmp.split("/")
						tree_name = parts[-1]		
						tree_name_id = tree_name.split(".nhx")[0]

						#tree.show()
						#gene tree build by Esaie's approach
						tree.write(outfile = "ensemblGeneTree/treeForComparaisonCorrectedEnsembl/" + tree_name_id + ".nw", format=9)

			
						for node in tree.traverse("postorder"):
							if node.is_leaf():				
								name = node.name
								node.name = transciptToGene[name]
								if node.name in ensembl_leaves:
									pass
								else:
									node.delete()					


						#tree.show()
						tree.write(outfile = "ensemblGeneTree/treeForComparaisonCorrectedEnsembl2/" + tree_name_id + ".nw", format=9)
			
						corrected_trees_by_ortho_group_nj.append(tree.write(format=9))

						leaves_Of_corrected_tree = tree.get_leaf_names()




					#Ensembl gene Trees
					ensembl_tree = Tree("ensemblGeneTree/treeForComparaisonEnsembl2/" + x + ".nw")
		

					while len(ensembl_tree.get_children()) == 1:
						ensembl_tree = ensembl_tree.get_children()[0]

					for node in ensembl_tree.traverse("preorder"):
						if node.is_leaf():
							if node.name in leaves_Of_corrected_tree:
								pass
							else:
								node.delete()
						else:
							if len(node.get_children()) == 1:												
								if not(node.is_root()):	
									parent = node.get_ancestors()[0]
									child = node.get_children()[0]
									child2 = deepcopy(child)
									child.delete()						
									parent.add_child(child2)
									node.delete()						
	
					while len(ensembl_tree.get_children()) == 1:
						ensembl_tree = ensembl_tree.get_children()[0]
			
					for node in ensembl_tree.traverse("preorder"):
						if node.is_leaf():
							if node.name in leaves_Of_corrected_tree:
								pass
							else:
								node.delete()
						else:
							if len(node.get_children()) == 1:												
								if not(node.is_root()):	
									parent = node.get_ancestors()[0]
									child = node.get_children()[0]
									child2 = deepcopy(child)
									child.delete()						
									parent.add_child(child2)
									node.delete()	

					print(ensembl_tree)

					ensembl_tree.write(outfile = "ensemblGeneTree/treeForComparaisonEnsembl2/" + x + ".nw", format=9)

					all_trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "w")

					for tree in corrected_trees_by_ortho_group:		
						all_trees.write(tree + "\n")

					for tree in corrected_trees_by_ortho_group_nj:		
						all_trees.write(tree + "\n")

					all_trees.write(ensembl_tree.write(format=9) + "\n")

					nb_of_trees = len(corrected_trees_by_ortho_group) + len(corrected_trees_by_ortho_group_nj) + 1 

					all_trees.close()

					mafft_cline1 = MafftCommandline(input="ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByOrthoGroup.fasta")
					mafft_cline2 = MafftCommandline(input="ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByEnsemblGeneTree.fasta")
					stdout1, stderr1 = mafft_cline1()
					stdout2, stderr2 = mafft_cline2()


		
					with open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByOrthoGroup.fasta", "w") as handle: 
						 handle.write(stdout1)

					with open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByEnsemblGeneTree.fasta", "w") as handle: 
						 handle.write(stdout2)

					main_fasta_to_phylip_one("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByOrthoGroup.fasta", "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy")
					main_fasta_to_phylip_one("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByEnsemblGeneTree.fasta", "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy")
		


					os.system("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					os.system("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.txt")

					os.system("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk")

					os.system("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy"+ "._seq_OrthoGroup_consel")
		
					#print("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -n 1 -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					#print("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.txt")

					#print("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk")

					#print("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy"+ "._seq_OrthoGroup_consel")

					#exit()
					"""
					print("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_2_Aln.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -n 1 -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					print("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_2_Aln.phy" + "_phyml_lk.txt")

					print("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_2_Aln.phy" + "_phyml_lk")

					print("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_2_Aln.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_2_Aln.phy"+ "._seq_Ensembl_consel")			
					"""
					"""
					os.system("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					os.system("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + "_phyml_lk.txt")

					os.system("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + "_phyml_lk")

					os.system("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy"+ "._seq_Ensembl_consel")			
					"""
					print(cmpt)
			except:
				pass

def main_corrected_reconstructEnsemblTree():	
	error_familly = open("ensemblGeneTree/error_familly.txt", "a")
	cmpt = 0
	listes = [[1,0]]
	#listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]	
	#random.shuffle(listes)
	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
			
		for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_" + sub_path + "/*_root.nhx"): 
			try:
				#for x in glob.glob("clusters_" + sub_path + "/fuzzyCMeans/ENSGT00390000016811_cluster_1_root.nhx"): 
				print(x)
				cmpt += 1

		
				mapping_familly = open("ensemblGeneTree/missing_familly.txt", "a")			
				x2 =x

				file = open(x, "r")
				x = str(x)
				x= x.split("/")
				x = x[-1]		
				x = x.split("_")[0]
				
				trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "r")
				lines = trees.readlines()
				trees.close()
				if len(lines)>3:
					pass
				else:	  			
					tree_build_by_ensembl = x2



					mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
					mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")

					if not (os.path.exists("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt")):
						mapping_familly.write(x + "\n")
						mapping_familly.close()
						continue
					mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")


					mapping_species = {}
					mapping_genes = {}
					geneToseq = {}
					transciptToGene = {}
					mapping_transcript_gene_protein_from_emsenbl = {}
					ensembl_leaves = []

					lines = mapping_transcript_gene_protein_from_emsenbl_file.readlines()
					for line in lines:
						line = line.replace("\n", "")
						parts = line.split("\t")
						if len(parts) == 4:
							mapping_transcript_gene_protein_from_emsenbl[parts[1]] = [parts[0], parts[3], parts[2]]
							ensembl_leaves.append(parts[1])


					lines = mapping_species_file.readlines()	

					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						#print(line)
						mapping_species[line[0]] = line[1]

					lines = mapping_genes_file.readlines()	
					#print(lines)
					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						mapping_genes[line[0]] = [line[1], line[2], line[3]]
						transciptToGene[line[1]] = line[2]


	

					file = open(tree_build_by_ensembl, "r")
					tree_str = file.read()
					tree_str = tree_str.replace("\n", "")
					tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
					#print(tree_str)
					trees_build_by_ensembl_nw = Tree(tree_str)
					for node in trees_build_by_ensembl_nw.traverse("postorder"):
						if node.is_leaf():				
							name = node.name
							node.name = node.name.split("_")[0]

					all_trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "a")


					all_trees.write(trees_build_by_ensembl_nw.write(format=9) + "\n")


					all_trees.close()

					#mafft_cline1 = MafftCommandline(input="ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByOrthoGroup.fasta")
					#mafft_cline2 = MafftCommandline(input="ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByEnsemblGeneTree.fasta")
					#stdout1, stderr1 = mafft_cline1()
					#stdout2, stderr2 = mafft_cline2()


		
					#with open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByOrthoGroup.fasta", "w") as handle: 
					#	 handle.write(stdout1)

					#with open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByEnsemblGeneTree.fasta", "w") as handle: 
					#	 handle.write(stdout2)
					#print("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByOrthoGroup.fasta", "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy")
					#main_fasta_to_phylip_one("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByOrthoGroup.fasta", "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy")
					#main_fasta_to_phylip_one("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByEnsemblGeneTree.fasta", "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy")
		


					os.system("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					os.system("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.txt")

					os.system("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk")

					os.system("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy"+ "._seq_OrthoGroup_consel")

					print(cmpt)	
			except:
				pass
		
		
###############################################################################
class Sequence(object):
	"""The Sequence object has a string *header* and
	various representations."""

	def __init__(self, header, seq):
		self.header = re.findall('^>(\S+)', header)[0]
		self.seq = seq

	def __len__(self):
		return len(self.seq)

	@property
	def phylip(self):
		return self.header + " " + self.seq.replace('.','-') + "\n"

	@property
	def fasta(self):
		return ">" + self.header + "\n" + self.seq + "\n"

def fasta_parse(path):
	"""Reads the file at *path* and yields
	   Sequence objects in a lazy fashion"""
	header = ''
	seq = ''
	with open(path) as f:
		for line in f:
			line = line.strip('\n')
			if line.startswith('>'):
				if header: yield Sequence(header, seq)
				header = line
				seq = ''
				continue
			seq += line
	yield Sequence(header, seq)


def main_fasta_to_phylip_one(fa_path, ph_path):
	# Check that the path is valid #
	if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
	# Use our two functions #
	seqs = fasta_parse(fa_path)
	# Write the output to temporary file #
	tm_path = ph_path + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
	# Count the sequences #
	count = 0
	with open(tm_path, 'w') as f:
		for seq in seqs:
			f.write(seq.phylip)
			count += 1
	# Add number of entries and length at the top #
	with open(tm_path, 'r') as old, open(ph_path, 'w') as new:
		new.write(" " + str(count) + " " + str(len(seq)) + "\n")
		new.writelines(old)
	# Clean up #
	command = "phyml -i "  + ph_path
	os.system(command)
		
	os.remove(tm_path)



def main_fasta_to_phylip_one(fa_path, ph_path):
	# Check that the path is valid #
	if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
	# Use our two functions #
	seqs = fasta_parse(fa_path)
	# Write the output to temporary file #
	tm_path = ph_path + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
	# Count the sequences #
	count = 0
	with open(tm_path, 'w') as f:
		for seq in seqs:
			f.write(seq.phylip)
			count += 1
	# Add number of entries and length at the top #
	with open(tm_path, 'r') as old, open(ph_path, 'w') as new:
		new.write(" " + str(count) + " " + str(len(seq)) + "\n")
		new.writelines(old)

		
	os.remove(tm_path)


def main_corrected_reconstructIsosel_withIsoselSeq():	
	error_familly = open("ensemblGeneTree/error_familly.txt", "a")
	cmpt = 0
	listes = [[1,0]]
	#listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]	
	#random.shuffle(listes)
	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
			
		for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_" + sub_path + "/*_root.nhx"): 
			#try:
				
			print(x)
			cmpt += 1
	
			mapping_familly = open("ensemblGeneTree/missing_familly.txt", "a")			
			x2 =x

			file = open(x, "r")
			x = str(x)
			x= x.split("/")
			x = x[-1]		
			x = x.split("_")[0]
			
			trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "r")
			lines = trees.readlines()
			trees.close()
			if len(lines) != 5:
				pass
			else:	  			
				transcrit2gene = {}
				transcrit2geneFile = open("initialSource/" + x + "_initialsource2target.txt", "r")
				lines = transcrit2geneFile.readlines()
				for l in lines:
					l = l.replace("\n", "")
					p = l.split(" ")
					cdsid = p[0]
					geneid = p[1]
					transcrit2gene[cdsid] = geneid
				
				isosel_fasta = open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByIsoSel.fasta", "w")					

				treef = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "r")
				lines = treef.readlines()
				l = lines[0]
				l = l.replace("\n", "")
				tree = Tree(l)
				gene_used = tree.get_leaf_names()				   
				
				initial_seq = {}
				
				for record in SeqIO.parse("initialSource/"+ x + "_initialsource.fasta", "fasta"):
					initial_seq[record.id] = record.seq

				for record in SeqIO.parse("isoselSeq/"+ x + "_filtered.fasta", "fasta"):
					cdsid = record.id

					if transcrit2gene[cdsid] in gene_used:
						seq = str(initial_seq[cdsid])
						isosel_fasta.write(">" + transcrit2gene[cdsid] + "\n")
						isosel_fasta.write(seq + "\n")

				isosel_fasta.close()
		  


				mafft_cline1 = MafftCommandline(input="ensemblGeneTree/sequences_" + sub_path + "/" + x + "_ByIsoSel.fasta")
				stdout1, stderr1 = mafft_cline1()



	
				with open("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByIsoSel.fasta", "w") as handle: 
					 handle.write(stdout1)


				main_fasta_to_phylip_one("ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Alni_ByIsoSel.fasta", "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByIsoSel.phy")

				os.system("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByIsoSel.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -o lr --print_site_lnl --no_memory_check --quiet -d aa")

				os.system("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByIsoSel.phy" + "_phyml_lk.txt")

				os.system("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByIsoSel.phy" + "_phyml_lk")

				os.system("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByIsoSel.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByIsoSel.phy"+ "._seq_ByIsoSel_consel")

				print(cmpt)	
			#except Exception as e:
				#print(e)
				#pass

def main_corrected_reconstruct_withEsemblSeq():	
	error_familly = open("ensemblGeneTree/error_familly.txt", "a")
	cmpt = 0
	listes = [[1,0]]
	#listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]	
	#random.shuffle(listes)
	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
			
		for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_" + sub_path + "/*_root.nhx"): 
			try:
				
				print(x)
				cmpt += 1
		
				mapping_familly = open("ensemblGeneTree/missing_familly.txt", "a")			
				x2 =x

				file = open(x, "r")
				x = str(x)
				x= x.split("/")
				x = x[-1]		
				x = x.split("_")[0]
				
				trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "r")
				lines = trees.readlines()
				trees.close()
				if len(lines) != 5:
					pass
				else:	  			
					tree_build_by_isosel = "isoselSeq/" + x + "_root.nhx"



					mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
					mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")

					if not (os.path.exists("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt")):
						mapping_familly.write(x + "\n")
						mapping_familly.close()
						continue
					mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")


					mapping_species = {}
					mapping_genes = {}
					mapping_genes2 = {}
					geneToseq = {}
					transciptToGene = {}
					mapping_transcript_gene_protein_from_emsenbl = {}
					ensembl_leaves = []

					lines = mapping_transcript_gene_protein_from_emsenbl_file.readlines()
					for line in lines:
						line = line.replace("\n", "")
						parts = line.split("\t")
						if len(parts) == 4:
							mapping_transcript_gene_protein_from_emsenbl[parts[1]] = [parts[0], parts[3], parts[2]]
							ensembl_leaves.append(parts[1])


					lines = mapping_species_file.readlines()	

					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						#print(line)
						mapping_species[line[0]] = line[1]

					lines = mapping_genes_file.readlines()	
					#print(lines)
					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						mapping_genes[line[0]] = [line[1], line[2], line[3]]
						transciptToGene[line[1]] = line[2]
						var = line[0]
						var = var.split("_")[0]
						mapping_genes2[var] = [line[1], line[2], line[3]]


	

					file = open(tree_build_by_isosel, "r")
					tree_str = file.read()
					tree_str = tree_str.replace("\n", "")
					tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
					#print(tree_str)
					tree_build_by_isosel_nw = Tree(tree_str)
					for node in tree_build_by_isosel_nw.traverse("postorder"):
						if node.is_leaf():				
							name = node.name
							node.name = mapping_genes2[node.name.split("_")[0]][1]
					print(tree_build_by_isosel_nw.write(format=9))
					#print(mapping_genes)
					#exit()
					#all_trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "a")


					#all_trees.write(tree_build_by_isosel_nw.write(format=9) + "\n")


					#all_trees.close()

			  


					os.system("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					os.system("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + "_phyml_lk.txt")

					os.system("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + "_phyml_lk")

					os.system("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByEnsemblGeneTree.phy"+ "._seq_ByEnsembl_consel")

					print(cmpt)	
			except Exception as e:
				print(e)
				pass
		
def main_corrected_reconstructIsosel():	
	error_familly = open("ensemblGeneTree/error_familly.txt", "a")
	cmpt = 0
	listes = [[1,0]]
	#listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]	
	#random.shuffle(listes)
	for e in  listes:
		weightStructure = e[0]
		weightSequence = e[1]
		sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
			
		for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_" + sub_path + "/*_root.nhx"): 
			try:
				
				print(x)
				cmpt += 1
		
				mapping_familly = open("ensemblGeneTree/missing_familly.txt", "a")			
				x2 =x

				file = open(x, "r")
				x = str(x)
				x= x.split("/")
				x = x[-1]		
				x = x.split("_")[0]
				
				trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "r")
				lines = trees.readlines()
				trees.close()
				if len(lines) != 4:
					pass
				else:	  			
					tree_build_by_isosel = "isoselSeq/" + x + "_root.nhx"



					mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
					mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")

					if not (os.path.exists("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt")):
						mapping_familly.write(x + "\n")
						mapping_familly.close()
						continue
					mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")


					mapping_species = {}
					mapping_genes = {}
					mapping_genes2 = {}
					geneToseq = {}
					transciptToGene = {}
					mapping_transcript_gene_protein_from_emsenbl = {}
					ensembl_leaves = []

					lines = mapping_transcript_gene_protein_from_emsenbl_file.readlines()
					for line in lines:
						line = line.replace("\n", "")
						parts = line.split("\t")
						if len(parts) == 4:
							mapping_transcript_gene_protein_from_emsenbl[parts[1]] = [parts[0], parts[3], parts[2]]
							ensembl_leaves.append(parts[1])


					lines = mapping_species_file.readlines()	

					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						#print(line)
						mapping_species[line[0]] = line[1]

					lines = mapping_genes_file.readlines()	
					#print(lines)
					for line in lines:
						line = line.replace("\n", "")
						line = line.split("\t")
						mapping_genes[line[0]] = [line[1], line[2], line[3]]
						transciptToGene[line[1]] = line[2]
						var = line[0]
						var = var.split("_")[0]
						mapping_genes2[var] = [line[1], line[2], line[3]]


	

					file = open(tree_build_by_isosel, "r")
					tree_str = file.read()
					tree_str = tree_str.replace("\n", "")
					tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
					#print(tree_str)
					tree_build_by_isosel_nw = Tree(tree_str)
					for node in tree_build_by_isosel_nw.traverse("postorder"):
						if node.is_leaf():				
							name = node.name
							node.name = mapping_genes2[node.name.split("_")[0]][1]
					print(tree_build_by_isosel_nw.write(format=9))
					#print(mapping_genes)
					#exit()
					all_trees = open("ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw", "a")


					all_trees.write(tree_build_by_isosel_nw.write(format=9) + "\n")


					all_trees.close()

			  


					os.system("PhyML-3.1/PhyML-3.1_linux64 -i " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + " -u " + "ensemblGeneTree/trees_" + sub_path + "/" + x + ".nw"  + " -o lr --print_site_lnl --no_memory_check --quiet -d aa")

					os.system("ensemblGeneTree/consel/bin/makermt --phyml " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.txt")

					os.system("ensemblGeneTree/consel/bin/consel "+ "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk")

					os.system("ensemblGeneTree/consel/bin/catpv " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy" + "_phyml_lk.pv -s 1 > " + "ensemblGeneTree/sequences_" + sub_path + "/" + x + "_Aln_ByOrthoGroup.phy"+ "._seq_OrthoGroup_consel")

					print(cmpt)	
			except:
				pass
		
###############################################################################
class Sequence(object):
	"""The Sequence object has a string *header* and
	various representations."""

	def __init__(self, header, seq):
		self.header = re.findall('^>(\S+)', header)[0]
		self.seq = seq

	def __len__(self):
		return len(self.seq)

	@property
	def phylip(self):
		return self.header + " " + self.seq.replace('.','-') + "\n"

	@property
	def fasta(self):
		return ">" + self.header + "\n" + self.seq + "\n"

def fasta_parse(path):
	"""Reads the file at *path* and yields
	   Sequence objects in a lazy fashion"""
	header = ''
	seq = ''
	with open(path) as f:
		for line in f:
			line = line.strip('\n')
			if line.startswith('>'):
				if header: yield Sequence(header, seq)
				header = line
				seq = ''
				continue
			seq += line
	yield Sequence(header, seq)


def main_fasta_to_phylip_one(fa_path, ph_path):
	# Check that the path is valid #
	if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
	# Use our two functions #
	seqs = fasta_parse(fa_path)
	# Write the output to temporary file #
	tm_path = ph_path + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
	# Count the sequences #
	count = 0
	with open(tm_path, 'w') as f:
		for seq in seqs:
			f.write(seq.phylip)
			count += 1
	# Add number of entries and length at the top #
	with open(tm_path, 'r') as old, open(ph_path, 'w') as new:
		new.write(" " + str(count) + " " + str(len(seq)) + "\n")
		new.writelines(old)
	# Clean up #
	command = "phyml -i "  + ph_path
	os.system(command)
		
	os.remove(tm_path)



def main_fasta_to_phylip_one(fa_path, ph_path):
	# Check that the path is valid #
	if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
	# Use our two functions #
	seqs = fasta_parse(fa_path)
	# Write the output to temporary file #
	tm_path = ph_path + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
	# Count the sequences #
	count = 0
	with open(tm_path, 'w') as f:
		for seq in seqs:
			f.write(seq.phylip)
			count += 1
	# Add number of entries and length at the top #
	with open(tm_path, 'r') as old, open(ph_path, 'w') as new:
		new.write(" " + str(count) + " " + str(len(seq)) + "\n")
		new.writelines(old)

		
	os.remove(tm_path)


def main_fasta_to_phylip(clusters_path):  
	for root, dirs, files in os.walk(clusters_path):  
		for filename in  files:
			command = "java -jar macse_v2.03.jar -prog alignSequences -seq "  + root + "/" + filename 
			os.system(command)		
			#print(filename)
			if filename.endswith(".fasta"):
				newname = filename.split(".")[0]
				main_fasta_to_phylip_one(clusters_path  +"/" +filename, clusters_path + "/phylipFormat/" +newname.split("_")[0]+".fa")



def main_stat():	
	error_familly = open("ensemblGeneTree/error_familly.txt", "a")
	cmpt = 0	
	all_stat = []
	for path in ["0_100", "20_80", "40_60", "50_50", "60_40", "80_20", "100_0"]:
		stats = []
		for x in glob.glob("clusters_" + path + "/fuzzyCMeans/*_nj_1_root.nhx"): 
		#for x in glob.glob("clusters_" + sub_path + "/fuzzyCMeans/ENSGT00390000016811_cluster_1_root.nhx"): 
			print(x)
			cmpt += 1			
			
			mapping_familly = open("ensemblGeneTree/missing_familly.txt", "a")			
			x2 =x

			file = open(x, "r")
			x = str(x)
			x= x.split("/")
			x = x[-1]		
			x = x.split("_")[0]
				

			trees_build_by_orthoGroup = []
			
			i = 1
			while os.path.exists("clusters_" + path + "/fuzzyCMeans/" + x + "_cluster_" + str(i) + "_root.nhx"):
				trees_build_by_orthoGroup.append("clusters_" + path + "/fuzzyCMeans/" + x + "_cluster_" + str(i) + "_root.nhx")
				break
				i += 1


			sequencesUsedByOrthoGroup = {}
			for record in SeqIO.parse("initialSource/" + x + "_initialsource.fasta", "fasta"):
				sequencesUsedByOrthoGroup[record.id] = str(record.seq)

			sequencesUsedByEnsemblGeneTree = {}
			for record in SeqIO.parse("ensemblGeneTree/seqFastaUsedByEnsembl/" + x + ".fasta", "fasta"):
				sequencesUsedByEnsemblGeneTree[record.id] = str(record.seq)

			source2target_dict = {}
			source2target = open("initialSource/" + x + "_initialsource2target.txt", "r")
			lines = source2target.readlines()
			for line in lines:
				line = line.replace("\n", "")
				parts = line.split(" ")
				id_transcript = parts[0]
				id_gene = parts[1]
				if id_gene in source2target_dict.keys():
					source2target_dict[id_gene].append(id_transcript)
				else:
					source2target_dict[id_gene] = [id_transcript]

			mapping_species_file = open("tmp/"+x+"_mapping_species.txt", "r")
			mapping_genes_file = open("tmp/"+x+"_mapping_transcript.txt", "r")
			if not (os.path.exists("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt")):
				mapping_familly.write(x + "\n")
				mapping_familly.close()
				continue
			mapping_transcript_gene_protein_from_emsenbl_file = open("ensemblGeneTree/transcript_protein_gene/" + x + "_update.txt", "r")


			mapping_species = {}
			mapping_genes = {}
			geneToseq = {}
			transciptToGene = {}
			mapping_transcript_gene_protein_from_emsenbl = {}
			ensembl_leaves = []
			
			lines = mapping_transcript_gene_protein_from_emsenbl_file.readlines()
			for line in lines:
				line = line.replace("\n", "")
				parts = line.split("\t")
				if len(parts) == 4:
					if parts[1] in source2target_dict.keys() and  len(source2target_dict[parts[1]])>1:
						mapping_transcript_gene_protein_from_emsenbl[parts[1]] = [parts[0], parts[3], parts[2]]
						ensembl_leaves.append(parts[1])


			lines = mapping_species_file.readlines()	

			for line in lines:
				line = line.replace("\n", "")
				line = line.split("\t")
				#print(line)
				mapping_species[line[0]] = line[1]

			lines = mapping_genes_file.readlines()	
			#print(lines)
			for line in lines:
				line = line.replace("\n", "")
				line = line.split("\t")
				mapping_genes[line[0]] = [line[1], line[2], line[3]]
				transciptToGene[line[1]] = line[2]


			trees_build_by_orthoGroup_nw = []



			for tree_st in trees_build_by_orthoGroup:
				file = open(tree_st, "r")
				tree_str = file.read()
				tree_str = tree_str.replace("\n", "")
				tree_str =  re.sub(r"\[([A-Za-z0-9_&:=$-]+)\]", r"",tree_str)
				#print(tree_str)
				trees_build_by_orthoGroup_nw.append(Tree(tree_str))

			#ensembl_tree = Tree("ensemblGeneTree/treeForComparaisonEnsembl2/" + x + ".nw")
			#ensembl_leaves = ensembl_tree.get_leaf_names()

			#Here I select just one tree to extract the set of used gene Ids.
			tree = trees_build_by_orthoGroup_nw[0]

			leavesNames = []
			selectedTranscriptUsedByOrthoGroup = []
			#print(tree.get_leaf_names())

			for node in tree.traverse("postorder"):
				if node.is_leaf():
					name = node.name
					name =  mapping_genes[name][1]
					#print(name)
					if name in ensembl_leaves and len(source2target_dict[name])>1:
						leavesNames.append(name)
						selectedTranscriptUsedByOrthoGroup.append(mapping_genes[node.name][0])

			#print(ensembl_leaves)
			#print(leavesNames)
			usedGeneIds = []
			for name in transciptToGene.keys():

					geneName = transciptToGene[name]
					
					if name in selectedTranscriptUsedByOrthoGroup:


						usedGeneIds.append(geneName)



			selectedTranscriptUsedByEnsemblGeneTree = []
			if geneName not in mapping_transcript_gene_protein_from_emsenbl.keys():
				continue				
			for geneName in usedGeneIds:
				proteinIdUsebForThisGeneID =  mapping_transcript_gene_protein_from_emsenbl[geneName][0]
				transcriptIdUsebForThisGeneID =  mapping_transcript_gene_protein_from_emsenbl[geneName][1]
				selectedTranscriptUsedByEnsemblGeneTree.append(transcriptIdUsebForThisGeneID)



			sharedTranscripts = [transcript for transcript in selectedTranscriptUsedByOrthoGroup if transcript in selectedTranscriptUsedByEnsemblGeneTree] 
			notSharedTranscripts1 = [transcript for transcript in selectedTranscriptUsedByOrthoGroup if transcript not in selectedTranscriptUsedByEnsemblGeneTree] 
			notSharedTranscripts2 = [transcript for transcript in selectedTranscriptUsedByEnsemblGeneTree if transcript not in selectedTranscriptUsedByOrthoGroup] 
			totalTranscripts = notSharedTranscripts1 + notSharedTranscripts2 + sharedTranscripts

			stat_file = open("ensemblGeneTree/sequences_" + path + "/" + x + "_stats_2.txt", "w")
			#stat_file.write(str((100.0*len(sharedTranscripts))/len(totalTranscripts)) + "\n")
			#stat_file.write("-".join(sharedTranscripts) +  "\t" + "-".join(notSharedTranscripts1) + "\t" + "-".join(notSharedTranscripts2) + "\t" + "-".join(totalTranscripts))
			stat_file.write(str(len(sharedTranscripts)) +  "\t"+ str(round((100.0*len(sharedTranscripts))/len(selectedTranscriptUsedByOrthoGroup),2)))
			stat_file.close()
			stats.append(round((100.0*len(sharedTranscripts))/len(selectedTranscriptUsedByOrthoGroup),2))

		
		all_stat.append(stats)
	c = "red"
	print(all_stat)
	print(len(all_stat))	
	plt.boxplot(all_stat,
						 vert=True,  # vertical box alignment
						 #patch_artist=True,  # fill with color
						 capprops=dict(color=c),
						 whiskerprops=dict(color=c),
						 flierprops=dict(color=c, markeredgecolor=c),
						 medianprops=dict(color=c),
						 labels=[r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks
	plt.xticks(fontsize=35)
	plt.yticks(fontsize=35)
	plt.ylabel('Average common transcripts', fontsize=35)
	plt.legend(prop={'size': 35})	
	plt.show()



#main_corrected_reconstructEnsemblTree()

try:
	#main_ensenmblTree()
	#main_corrected_ensemblTree()
	#main_corrected_reconstructEnsemblTree()
	main_corrected_reconstructIsosel()
except Exception as e:
	exc_type, exc_obj, exc_tb = sys.exc_info()
	fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
	print(exc_type, fname, exc_tb.tb_lineno)
	print(traceback.format_exc())

#main_stat()
#main_corrected_reconstruct_withEsemblSeq()
main_corrected_reconstructIsosel_withIsoselSeq()
"""
if __name__ == "__main__": 
	try:
		i =2
		#main_ensenmblTree()
		main_corrected_ensemblTree()
		#main_fasta_to_phylip("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/target")
	except :#Exception as e:
		#exc_type, exc_obj, exc_tb = sys.exc_info()
		#fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		#print(exc_type, fname, exc_tb.tb_lineno)
		#print(traceback.format_exc())
		pass
		
"""
