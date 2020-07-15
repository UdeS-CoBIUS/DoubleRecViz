from fusionMatrix import *
import glob
import numpy as np
import math
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
def correlation():
	correlation = []

	try:
		X1 = []
		Y1 = []
		X2 = []
		Y2 = []
		for structSim in glob.glob("similarityScores/microalignment_*_score.csv"): 
			ch = str(structSim)
			parts= ch.split("/")
			id = parts[-1].split("_")[1]
			#if id != "ENSGT00940000158566":
			#	continue
			seqSim = "similaritySeqScores/" +  id + "_microalignment.matrix.csv"
			#print(id)
			x_abs, y_ord, x_abs2, y_ord2 = correlation_matrix(seqSim, structSim)
			if len(x_abs) == 0:
				pass
			else:
				for x in x_abs:
					X1.append(x)
				for x in y_ord:
					Y1.append(x)
				for x in x_abs2:
					X2.append(x)
				for x in y_ord2:
					Y2.append(x)
			"""
			if cor != 100 and str(cor) != "nan":
				correlation.append(cor)
				print (id + "\t" + str(cor))		
			"""

		plt.scatter(X1, Y1, s=50, cmap='gray', c='black')	
		plt.scatter(X2, Y2, s=50, cmap='gray', c='red')
		#plt.axis('off')
		plt.show()
	except Exception as e:
		print(e)
	print("-----------------")
	print(np.mean(correlation))	
	print(max(correlation))
	print(min(correlation))
    

def correlation2():
	correlation = []

	try:
		X1 = []
		Y1 = []
		X2 = []
		Y2 = []
		for structSim in glob.glob("similarityScores/microalignment_*_score.csv"): 
			ch = str(structSim)
			parts= ch.split("/")
			id = parts[-1].split("_")[1]
			#if id != "ENSGT00940000158566":
			#	continue
			seqSim = "similaritySeqScores/" +  id + "_microalignment.matrix.csv"
			#print(id)
			ratio = correlation_matrix2(seqSim, structSim)
			correlation.append(ratio)
			print(np.mean(ratio))
	except Exception as e:
		print(e)
	print("-----------------")
	print(np.mean(correlation))	
	print(max(correlation))
	print(min(correlation))


def correlation3():
	correlation = []

	try:

		for structSim in glob.glob("similarityScores/microalignment_*_score.csv"):
			X1 = []
			Y1 = [] 
			ch = str(structSim)
			parts= ch.split("/")
			id = parts[-1].split("_")[1]
			#if id != "ENSGT00940000158566":
			#	continue
			seqSim = "similaritySeqScores/" +  id + "_microalignment.matrix.csv"
			#print(id)
			X, Y = correlation_matrix3(seqSim, structSim)
			cor = pearsonr(X,Y)
			correlation.append(cor[0])
            
			for x in X:
				X1.append(x)
			for y in Y:
				Y1.append(y)
			"""
			if len(X)>9:				
				#print(X1)			
				#print(Y1)
				tips = pd.DataFrame(data={"X": list(X1), "Y":list(Y1)	})
				#print(tips)
				sns.set(style="darkgrid")

				#tips = sns.load_dataset("tips")
				g = sns.jointplot("X", "Y", data=tips, kind="reg", xlim=(0,1), ylim=(0, 1), color="m", height=7)		
				g.set_axis_labels("sequence similarity", "splicing structure similarity")

				# Set title
				g.fig.suptitle("sequence similarity vs splicing structure similarity")

				# Format nicely.
				g.fig.tight_layout()

				# Reduce plot to make room for suptitle
				g.fig.subplots_adjust(top=0.90)

				# Save and close the figure
				plt.savefig("tes.png")

				#plt.close('all') 		
				#plt.scatter(X1, Y1)
				plt.show()
			"""
	except Exception as e:
		print(e)
	print(correlation)
	#print(np.mean(correlation))	


correlation3()
