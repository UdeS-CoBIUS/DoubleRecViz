import os
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob
import numpy as np
import matplotlib.pyplot as plt

nb_total_files = 0.0
for root, dirs, files in os.walk("ensemblGeneTree/sequences_40_60"): 
        dic = {}		
        for file in files:
            tmp = file.split("_")
            if file.count("_")>0 and len(tmp) > 0 :
                nb_total_files += 1        
                geneFamilyId = tmp[0]			
                prot_gene_spe = open("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/transcript_protein_gene/" + geneFamilyId + "_update.txt", "r")
                lines =  prot_gene_spe.readlines()
                nb_total_lines = 0.0
                nb_includes_lines = 0.0
                for l in lines:
                    parts = l.split("\t")
                    if len(parts) == 4:
                        nb_includes_lines += 1
                    nb_total_lines += 1
    
                ratio = nb_includes_lines/nb_total_lines
                ratio = int(ratio*10)
                if ratio == 10:
                    ratio = 9
                if ratio in dic.keys():
                    dic[ratio] += 1
                else:
                    dic[ratio] = 1
                #print(dic)  
                
for k, v in  dic.items():
    dic[k] = v*100/nb_total_files

print(dic)    

barWidth = 0.6

# Choose the height of the blue bars 
#moyenne equal
bars1 = []
for i in range(10):
    bars1.append(dic[i])

# The x position of bars
#r1 = [x + 0.1 for x in np.arange(len(bars1))]
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + 0.2+barWidth for x in r1]

# Create blue bars
rects1 = plt.bar(r2, bars1, width = barWidth, color = '#4285F4', edgecolor = '#4285F4', capsize=40)

# general layout
plt.xticks([r + barWidth for r in range(len(bars1))], ["[0, 0.1[", "[0.1, 0.2[", "[0.2, 0.3[", "[0.3, 0.4[", "[0.4, 0.5[", "[0.5, 0.6[", "[0.6, 0.7[", "[0.7, 0.8[", "[0.8, 0.9[", "[0.9, 1]"])
plt.xticks(fontsize=60, rotation=50)
plt.yticks(fontsize=60)
plt.legend(prop={'size': 75})
plt.xlabel('Ratio value classes', fontsize=90)
plt.ylabel('Ratio percentage', fontsize=70)
#plt.savefig("RI.png")
# Show graphic
def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()        
        plt.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                 str(height)+  "%",
                ha='center', va='bottom',  fontsize=30)
 
#autolabel(rects1)                
#autolabel(rects2)                
#autolabel(rects3)                
plt.show()

