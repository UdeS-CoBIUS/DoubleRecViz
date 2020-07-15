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
import matplotlib.pyplot as plt
import pandas as pd

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
			try:
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
			except:
			    pass
		
		all_stat.append(stats)
	c = "black"
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
	plt.xticks(fontsize=50, rotation=50)
	plt.yticks(fontsize=50)
	plt.ylabel('Average cxommon transcripts', fontsize=50)
	plt.legend(prop={'size': 60})	
	plt.show()

main_stat()
