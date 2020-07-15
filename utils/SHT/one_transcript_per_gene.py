import glob
import os
dict = {}
for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/initialSource/*initialsource2target.txt"): 
    tmp = x.split("/")
    tmp = tmp[-1]  
    tmp = tmp.split("_")
    id = tmp[0]
    if os.path.exists("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/clusters_0_100/fuzzyCMeans/" + id + "_cluster_1_root.nhx"):
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
        nb_gene_tmp = 0
        for k, v in tmp_dict.items():
            if v == 1:
                nb_gene_tmp += 1
        fraction_gene_with_one_transcript = nb_gene_tmp*100.0/len(tmp_dict.keys())

        if fraction_gene_with_one_transcript in dict.keys():
            dict[fraction_gene_with_one_transcript] += 1
        else:
            dict[fraction_gene_with_one_transcript] = 1                    

print(dict)
dic_f = {}
total = 0.0
for i, v in dict.items():  
    i = int(i)
    total += v
    if i in dic_f.keys():  
        dic_f[i] += v
    else:
        dic_f[i] = v

for i, v in dic_f.items():  
    dic_f[i] = 100*dic_f[i]/total

dic_f2 = {}
for i, v in dic_f.items():  
    i = int(i/10)
    total += v
    if i in dic_f2.keys():  
        dic_f2[i] += v
    else:
        dic_f2[i] = v    
print(dic_f2)

dic_f2[1] += dic_f2[0]
del dic_f2[0]

keys = sorted(dic_f2.keys(), reverse=True)
bars1 = []
for k in keys:
    bars1.append(dic_f2[k])
print(bars1)   


         

# libraries
import numpy as np
import matplotlib.pyplot as plt

# width of the bars
barWidth = 0.6

# Choose the height of the blue bars 
#moyenne equal
#bars1 = [4.7244094488,15.7480314961,14.1732283465,23.6220472441,14.9606299213,6.2992125984,11.0236220472,3.937007874,5.5118110236]



# The x position of bars
#r1 = [x + 0.1 for x in np.arange(len(bars1))]
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + 0.2+barWidth for x in r1]

# Create blue bars
rects1 = plt.bar(r2, bars1, width = barWidth, color = '#4285F4', edgecolor = '#4285F4', capsize=40)


# general layout
plt.xticks([r + barWidth for r in range(len(bars1))], ["[100, 90[", "[90, 80[", "[80, 70[", "[70, 60[", "[60, 50[", "[50, 40[", "[40, 30[", "[30, 20[", "[20, 10[", "[10, 0]"])
plt.xticks(fontsize=60, rotation=50)
plt.yticks(fontsize=60)
plt.legend(prop={'size': 40})
plt.xlabel('Percentage of gene families per percentage \n of genes having a single CDS', fontsize=70)
plt.ylabel('Percentage of gene families', fontsize=40)

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()        
        plt.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                 str(height)+  "%",
                ha='center', va='bottom',  fontsize=30)
             
plt.show()

