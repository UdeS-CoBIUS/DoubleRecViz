import glob
import os
dict = {}
values = ["0_100", "20_80", "40_60", "50_50", "60_40", "80_20", "100_0"]
for path in values:	
	nb_isoform_per_gene = {}
	nb_gene_tmp = 0
	for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/initialSource/*initialsource2target.txt"): 
		tmp = x.split("/")
		tmp = tmp[-1]  
		tmp = tmp.split("_")
		id = tmp[0]
		keep = False
		if os.path.exists("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/clusters_" + path + "/fuzzyCMeans/" + id + "_cluster_1_root.nhx"):
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
					keep = True
				else:
					tmp_dict[gene_id] = 1
			

			if keep:
				for k, v in tmp_dict.items():
					nb_gene_tmp += 1
					if v in nb_isoform_per_gene.keys():
						nb_isoform_per_gene[v] +=1
					else:
						nb_isoform_per_gene[v] =1

	for k, v in nb_isoform_per_gene.items():
		nb_isoform_per_gene[k] = v*100/nb_gene_tmp

	for i in sorted (nb_isoform_per_gene.keys()) :  
		print(i, nb_isoform_per_gene[i])	   
			 
