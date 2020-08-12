# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

'''

import sys
import os
import os.path
from prune_species_tree import *
from copy_gene_tree_init import *
from compute_gene_trees import *
from compute_transcript_trees import *
from reconcile_trees import *
from NHXtoDoubleRecPhyloXML import *
from functions import *

Dir = sys.argv[1]


print("Generating test data from "+Dir+" ...")
print("")

if Dir[-1] == "/":
    Dir = Dir[0:-1]

# gene famiky name
f = Dir.split("/")[-1]

Dir += "/"
Transcript = Dir+"Sequences/"+f+"_transcript.fasta"
Gene = Dir+"Sequences/"+f+"_gene.fasta"
Transcript2gene = Dir+"Sequences/"+f+"_transcript2gene.txt"
Gene2species = Dir+"Sequences/"+f+"_gene2species.txt"
Genetree_init = Dir+"Sequences/"+f+"_genetree.nw"
Speciestree_init = Dir+"Sequences/"+f+"_speciestree.nw"


#****************************************************************
print("1. Computing trees")

Dir_trees = Dir + "Trees/"

os.system("mkdir " + Dir_trees)


#==========================================
print("... Species tree")
s = prune_species_tree(Speciestree_init,Gene2species)

pruned_speciestree = Dir_trees +f+"_pruned_speciestree.nw"
write_tree(s,pruned_speciestree,'w')
print("...... species tree written in :\n          ",pruned_speciestree)
print("")

speciestree = [pruned_speciestree]
#==========================================
print("... Gene trees")
if(os.path.exists(Genetree_init)):
    s = copy_gene_tree_init(Genetree_init)

    gene_tree_init = Dir_trees +f+"_genetree_init.nw"
    write_tree(s,gene_tree_init,'w')
    print("...... init gene tree written in :\n          ",gene_tree_init)
    genetrees = [gene_tree_init]
else:
    genetrees = []
#==========================================
s_treebest, s_phyml = compute_gene_trees(Transcript,
                       Transcript2gene,
                       Gene2species,
                       pruned_speciestree)

gene_tree_treebest = Dir_trees +f+"_genetree_treebest.nw"
write_tree(s_treebest,gene_tree_treebest,'w')
print("...... treebest gene tree written in :\n          ",gene_tree_treebest)
genetrees.append(gene_tree_treebest)

gene_tree_phyml = Dir_trees +f+"_genetree_phyml.nw"
write_tree(s_phyml,gene_tree_phyml,'w')
print("...... phyml gene tree written in :\n          ",gene_tree_phyml)
genetrees.append(gene_tree_phyml)
print("")

#==========================================
print("... Transcript trees")
s_treebest_best_treebest, s_treebest_best_phyml, s_phyml = compute_transcript_trees(Transcript,Transcript2gene,gene_tree_treebest, gene_tree_phyml)

transcript_tree_treebest_treebest = Dir_trees +f+"_transcripttree_treebest_treebest.nw"
write_tree(s_treebest_best_treebest,transcript_tree_treebest_treebest,'w')
print("...... treebest transcript tree (using treebest gene tree) written in :\n          ",transcript_tree_treebest_treebest)
transcripttrees = [transcript_tree_treebest_treebest]

transcript_tree_treebest_phyml = Dir_trees +f+"_transcripttree_treebest_phyml.nw"
write_tree(s_treebest_best_phyml,transcript_tree_treebest_phyml,'w')
print("...... treebest transcript tree (using phyml gene tree) written in :\n          ",transcript_tree_treebest_phyml)
transcripttrees.append(transcript_tree_treebest_phyml)

transcript_tree_phyml = Dir_trees +f+"_transcripttree_phyml.nw"
write_tree(s_phyml,transcript_tree_phyml,'w')
print("...... phyml transcript tree written in :\n          ",transcript_tree_phyml)
transcripttrees.append(transcript_tree_phyml)
print("")

#****************************************************************
print("2. Computing Reconciliations")

Dir_reconciliations = Dir + "Reconciliations/"

os.system("mkdir " + Dir_reconciliations)

all_rec_NHX = Dir_reconciliations +f+"_all_rec.nhx"
os.system("rm "+all_rec_NHX)

i = 0
for gene_tree in genetrees:
    genespecies_rec = Dir_reconciliations +f+"_gene"+str(i)+"_species.nhx"
    s = reconcile(gene_tree,pruned_speciestree,Gene2species)
    write_tree(s,genespecies_rec,'w')
    write_tree(s,all_rec_NHX,'a')
    j = 0
    for transcript_tree in transcripttrees:
        transcriptgene_rec = Dir_reconciliations +f+"_transcript"+str(j)+"_gene"+str(i)+".nhx"
        s = reconcile(transcript_tree, gene_tree,Transcript2gene)
        write_tree(s,transcriptgene_rec,'w')
        write_tree(s,all_rec_NHX,'a')
        j += 1
    write_tree("//",all_rec_NHX,'a')
    i += 1
        
print("... Gene2species and Transcript2species reconciliations (nhx format) in :\n          ",Dir_reconciliations)
print("... All reconciliations (nhx format) in :\n          ",all_rec_NHX)
print("")

#****************************************************************
print("3. Convert from multiple NHX to DoubleRecPhyloXML")

#==========================================
print("... Converting from NHX to DoubleRecPhyloXML")

all_rec_XML = all_rec_NHX.split(".nhx")[0]+".xml"
s_xml = NHXtoDoubleRecPhyloXML(all_rec_NHX, pruned_speciestree)
write_tree(s_xml,all_rec_XML, 'w')
print("... All reconciliations (doubleRecPhyloXML format) in :\n          ",all_rec_XML)
print("")

os.system("rm filtalign.fa")
os.system("rm -rf __pycache__")
