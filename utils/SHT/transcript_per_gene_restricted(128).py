import os
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob
import numpy as np
import matplotlib.pyplot as plt

dict = {}
for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/initialSource/*initialsource2target.txt"): 
    flag = False
    file = open(x)
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
        for k, v in tmp_dict.items():
            if v in dict.keys():
                dict[v] += 1
            else:
                dict[v] = 1                    

nb_gene = 0
for k, v in dict.items():
    nb_gene += v           
X = []
Y = []    
for i in sorted (dict.keys()) :  
    print(i, dict[i], dict[i]*100/nb_gene)       
    X.append(i)
    Y.append(dict[i]*100/nb_gene)
print(X)    
print(Y)
         

barWidth = 0.6

# Choose the height of the blue bars 
#moyenne equal
bars1 = Y

# The x position of bars
#r1 = [x + 0.1 for x in np.arange(len(bars1))]
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + 0.2+barWidth for x in r1]

# Create blue bars
rects1 = plt.bar(r2, bars1, width = barWidth, color = '#4285F4', edgecolor = '#4285F4', capsize=40)

# general layout
plt.xticks([r + barWidth for r in range(len(bars1))], X)
plt.xticks(fontsize=60)
plt.yticks(fontsize=60)
plt.legend(prop={'size': 60})
plt.xlabel('Number of CDS', fontsize=90)
plt.ylabel('Percentage of genes', fontsize=80)
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

