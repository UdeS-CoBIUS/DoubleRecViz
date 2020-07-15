import os
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob
# libraries
import numpy as np
import matplotlib.pyplot as plt
#share_leaves_Boxplot.pdf
def sharedLeaves():
    all_list = []
    
    listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]

    for e in  listes:
        weightStructure = e[0]
        weightSequence = e[1]
        sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))
	    
        percentageSharedLeaves = []
        for x in glob.glob("ensemblGeneTree/sequences_" + sub_path + "/*_stats.txt"): 
            tmp = x.split("/")[-1]
            tmp = tmp.split("_")[0]
            intial_source2target = "initialSource/"+ tmp + "_initialsource2target.txt"
            flag = False
            file = open(intial_source2target, "r")
            lines = file.readlines()
            tmp_dict = {}
            for line in lines:
                line = line.replace("\n", "")
                parts = line.split(" ")
                gene_id = parts[1]
                cds_id = parts[0]        
                if gene_id in tmp_dict.keys():
                    tmp_dict[gene_id] += 1
                else:
                    tmp_dict[gene_id] = 1
            for k, v in tmp_dict.items():
                if v!=1:
                    flag = True
    
            if flag == False:
                pass
            else:           
                sharedLeavesFile = open(x, "r")
                lines = sharedLeavesFile.readlines()
                if len(lines)>0:
                    line = lines[0]
                    line = line.replace("\n", "")
                    parts = line.split("\t")
                    if len(parts)>1:
                        percentageSharedLeaves.append(float(parts[1]))

        all_list.append(percentageSharedLeaves)
                    

    c = "black"
    plt.boxplot(all_list,
	                     vert=True,  # vertical box alignment
	                     #patch_artist=True,  # fill with color
                         capprops=dict(color=c),
			             whiskerprops=dict(color=c),
			             flierprops=dict(color=c, markeredgecolor=c),
			             medianprops=dict(color=c),
	                     labels=[r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks
    plt.xticks(fontsize=70, rotation=50)
    plt.yticks(fontsize=70)
    plt.ylabel('Average common transcripts', fontsize=50)
    plt.legend(prop={'size': 60})	
    plt.show()

#computeStats()
sharedLeaves()
