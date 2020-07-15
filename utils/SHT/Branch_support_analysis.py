#!/usr/bin/python3.7
#		 Computing statistics on species tree branch supports
#
#		 Da Rocha Coimbra, Nilson Antonio
#		 nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#		 Ouangraoua, Aida
#		 aida.ouangraoua@usherbrooke.ca
#
#		 2020
#

#from Utils import *
from ete3 import Tree
import decimal
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

###############################################################################

# Directory for statistics files

###############################################################################

sub_path = "100_0"
species_tree_LT = "Astral/test_data/"+ sub_path + "_SpeciesTreeEnsembl.tre"
species_tree_LTBD = "Astral/test_data/"+ sub_path + "_SpeciesTreeReconstructEnsembl.tre"
species_tree_IsoSel = "Astral/test_data/"+ sub_path + "_SpeciesTreeIsoSel.tre"
species_tree_consensus_tree = "Astral/test_data/"+ sub_path + "_Species_consensus.tre"

files = [species_tree_consensus_tree, species_tree_LTBD, species_tree_LT, species_tree_IsoSel]

listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
#listes = [[0.5,0.5]]		
support_values_per_method = []

for e in  listes:
	weightStructure = e[0]
	weightSequence = e[1]
	sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))

	species_tree_SHT = "Astral/test_data/"+ sub_path + "_SpeciesTreeOrthoGroup.tre"
	
	files.append(species_tree_SHT)

print(files)

#files = [species_tree_consensus_tree, species_tree_LTBD, species_tree_IsoSel, species_tree_LT]
for f in files:
	support_values = []
	p = f.split("/")[-1]
	t_astral = Tree(f)
	outpufile = open("Astral/test_data/support_" + p, "w")
	for node in t_astral.traverse("postorder"):
		if( not node.is_leaf()):
			#nb += 1
			leaves = []
			for leaf in node:
				leaves.append(leaf.name)
			outpufile.write(str(node.support)+"\n")
			support_values.append(node.support)
	support_values_per_method.append(support_values)		
	outpufile.close()

print(files)
print(support_values_per_method)

#exit()	

###############################################################################

# Generate statistics on common clades between trees

###############################################################################
print("\nGenerating statistics on common clades between trees...")
start = time.time()
#files = [raxmlfinal,astrid20,astrid40,astrid50,astrid60,astrid80,astral20,astral40,astral50,astral60,astral80,astridfinal,astralfinal,overallfinal]

common_clades_stat_file = "Astral/test_data/common_clades.stats"
stat = open(common_clades_stat_file,"w")
line = ""
for file in files:
	method = file.split("/")[-1].split(".txt")[0]
	line += "\t"+method
stat.write(line+"\n")
for i in range(len(files)):
	filei = files[i]	
	ti = Tree(filei)
	method = filei.split("/")[-1].split(".txt")[0]
	line = method
	for j in range(i):
		line += "\t"
	line += "\t"
	for j in range(i):
		filej = files[j]
		file_open = open(filej, "r")
		content = file_open.readlines()[0]
		file_open.close()

		file_open = open(filej, "w")
		content = content.replace(":75.00:", ":", 10000)
		content = content.replace(":100.00:", ":", 10000)
		content = content.replace(":25.00:", ":", 10000)
		content = content.replace(":50.00:", ":", 10000)		
		file_open.write(content)
		file_open.close()		
		#filej = filej.replace(":75.00:", ":", 10000)
		print(filej, j)
		tj = Tree(filej, quoted_node_names=True, format=1)
		nb_ti = 0
		nb_tj = 0
		nb_common_clade = 0
		for node in ti.traverse("postorder"):
			if(not node.is_leaf()):
				nb_ti += 1
				leaves = []
				for leaf in node:
					leaves.append(leaf.name)
				if(tj.check_monophyly(values=leaves, target_attr="name", ignore_missing=True)[0]):
					nb_common_clade +=1
		for node in tj.traverse("postorder"):
			if( not node.is_leaf()):
				nb_tj += 1
		percent_common_clade = 100.0*nb_common_clade/min(nb_ti,nb_tj)
		line += "\t"+str(round(percent_common_clade,2))
	stat.write(line+"\n")
stat.close()
print("Statistics on common clades generated in "+str(time.time()-start)+" seconds")
print(common_clades_stat_file)

exit()
###############################################################################

# Generate statistics on branch supports

###############################################################################
print("\nGenerating statistics on branch supports...")

for i in range(len(files)):
	filei = files[i]  
	file_open = open(filei, "r")
	content = file_open.readlines()[0]
	file_open.close()
	print(filei)
	file_open = open(filei, "w")
	content = content.replace(":75.00:", ":", 10000)
	content = content.replace(":100.00:", ":", 10000)
	content = content.replace(":25.00:", ":", 10000)
	content = content.replace(":50.00:", ":", 10000)		
	file_open.write(content)
	file_open.close()		 
	ti = Tree(filei, quoted_node_names=True, format=1)

	file_name = filei.split("/")[-1]
	file_name = file_name.split(".")[0]

	support_overall_stat_file = "Astral/test_data/support_overall_"+file_name+".csv"
	outpufile = open(support_overall_stat_file, "w")
	outpufile.write("tree" + ","+ "support" +"," + "overall"+"\n")

	# 1
	nb = 0
	for node in ti.traverse("postorder"):
		if( not node.is_leaf()):
			nb += 1
			leaves = []
			for leaf in node:
				leaves.append(leaf.name)
			outpufile.write("astrid" + ","+ str(node.support)+"\n")
	outpufile.close()
	print("Statistics on branch support generated in ")
	


"""

###############################################################################

# Generate statistics on node repartition

###############################################################################
print("\nGenerating statistics on node repartition...")
start = time.time()

node_repartition_stat_file = statisticspath+ "/node_repartition.csv"
outpufile = open(node_repartition_stat_file, "w")
outpufile.write("tree" + ","+ "support"+"\n")
nb = 0
for node in t_overall.traverse("postorder"):
	if( not node.is_leaf()):
		nb += 1
		leaves = []
		for leaf in node:
			leaves.append(leaf.name)
		if(t_astrid.check_monophyly(values=leaves, target_attr="name")[0] and t_astral.check_monophyly(values=leaves, target_attr="name")[0] and t_raxml.check_monophyly(values=leaves, target_attr="name")[0]):
			outpufile.write("astrid-astral-raxml" + ","+ str(node.support)+"\n")
		elif(t_astrid.check_monophyly(values=leaves, target_attr="name")[0] and t_astral.check_monophyly(values=leaves, target_attr="name")[0]):
			outpufile.write("astrid-astral" + ","+ str(node.support)+"\n")
		elif(t_raxml.check_monophyly(values=leaves, target_attr="name")[0] and t_astral.check_monophyly(values=leaves, target_attr="name")[0]):
			outpufile.write("astral-raxml" + ","+ str(node.support)+"\n")
		elif(t_astrid.check_monophyly(values=leaves, target_attr="name")[0] and t_raxml.check_monophyly(values=leaves, target_attr="name")[0]):
			outpufile.write("astrid-raxml" + ","+ str(node.support)+"\n")
		else:
			outpufile.write("overall" + ","+ str(node.support)+"\n")
#print("overall",nb)
outpufile.close()
print("Statistics on node repartition generated in "+str(time.time()-start)+" seconds")
print(node_repartition_stat_file)

###############################################################################

# Generate branch support figure

###############################################################################
print("\nGenerating branch support figure...")
start = time.time()
support_overall = pd.read_csv(support_overall_stat_file)
node_repartition = pd.read_csv(node_repartition_stat_file)
support_figure_file = statisticspath+ "/support.png"

fig, axes = plt.subplots(2, 2,figsize=(12,12))
sns.boxplot(x="tree", y="support", hue = "overall", data=support_overall, ax = axes[0,0])
sns.countplot(x="tree", hue = "overall", data=support_overall, ax = axes[1,0])
sns.boxplot(x="tree", y="support", data=node_repartition, ax = axes[0,1])
sns.countplot(x="tree", data=node_repartition, ax = axes[1,1])
fig.savefig(support_figure_file)

print("Branch support figure generated in "+str(time.time()-start)+" seconds")
print(support_figure_file)
"""