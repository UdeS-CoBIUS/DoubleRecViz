from io import StringIO
from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
from itertools import permutations

def read_newick(treedata):
    handle = StringIO(treedata)
    return Phylo.read(handle, "newick")


listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]		
#listes = [[0.5,0.5]]		


for e in  listes:
	weightStructure = e[0]
	weightSequence = e[1]
	sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
	speciesTreeIsoSel_f = open("Astral/test_data/" + sub_path + "_SpeciesTreeIsoSel.tre", "r")
	speciesTreeReconstructEnsembl_f = open("Astral/test_data/" + sub_path + "_SpeciesTreeReconstructEnsembl.tre", "r")
	speciesTreeEnsembl_f = open("Astral/test_data/" + sub_path + "_SpeciesTreeEnsembl.tre", "r")
	speciesTreeOrthoGroup_f = open("Astral/test_data/" + sub_path + "_SpeciesTreeOrthoGroup.tre", "r")


	speciesTreeIsoSel_nw = Phylo.read("Astral/test_data/" + sub_path + "_SpeciesTreeIsoSel.tre", "newick")
	speciesTreeReconstructEnsembl_nw = Phylo.read("Astral/test_data/" + sub_path + "_SpeciesTreeReconstructEnsembl.tre", "newick") 
	speciesTreeEnsembl_nw = Phylo.read("Astral/test_data/" + sub_path + "_SpeciesTreeEnsembl.tre", "newick") 
	speciesTreeOrthoGroup_nw = Phylo.read("Astral/test_data/" + sub_path + "_SpeciesTreeOrthoGroup.tre", "newick") 


	trees = [speciesTreeIsoSel_nw, speciesTreeReconstructEnsembl_nw, speciesTreeEnsembl_nw, speciesTreeOrthoGroup_nw]
	majority_tree = majority_consensus(trees, 0.5)
	print('majority consensus for order:')
	#Phylo.draw_ascii(majority_tree)
	Phylo.write(majority_tree, "Astral/test_data/"+sub_path+"_Species_consensus.tre", 'newick')



"""
files = [species_tree_consensus_tree, species_tree_LTBD, species_tree_IsoSel, species_tree_LT]
for f in files:
	p = f.split("/")[-1]
	t_astral = Tree(f)
	outpufile = open("Astral/test_data/support_" + p, "w")
	for node in t_astral.traverse("postorder"):
		if( not node.is_leaf()):
			nb += 1
			leaves = []
			for leaf in node:
				leaves.append(leaf.name)
			outpufile.write("astral" + ","+ str(node.support)+"\n")

	outpufile.close()
"""