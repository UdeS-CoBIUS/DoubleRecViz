import glob
import os
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
def compute_alignment_muscle(nt_file, aa_file, nt_alignment, aa_alignment):
	print(nt_file, aa_file, nt_alignment, aa_alignment)
	
	command = "./muscle3.8.31_i86linux32 -in "  + aa_file + " -out " +  aa_alignment
	os.system(command)	

	nt_dict = {}
	for record in SeqIO.parse(nt_file, "fasta"):
		nt_dict[record.id] = str(record.seq)

	aa_aln_dict = {}
	for record in SeqIO.parse(aa_alignment, "fasta"):
		aa_aln_dict[record.id] = str(record.seq)		
	
	file = open(nt_alignment, "w")
	for id in aa_aln_dict.keys():
		if id in nt_dict.keys():
			back_translate_seq = ""
			aa_aln_seq = aa_aln_dict[id]
			nt_seq = nt_dict[id]
			nb_gap = 0
			for pos in range(len(aa_aln_seq)):
				if aa_aln_seq[pos] == "-":
					back_translate_seq += "---"
					nb_gap += 1
				else:
					pos_to_use = pos - nb_gap
					back_translate_seq += nt_seq[3*pos_to_use:3*pos_to_use+3]
			file.write(">" + id + "\n")
			file.write(back_translate_seq + "\n")

	file.close()


def recomputeEnsemblTree():
	for x in glob.glob("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_100_0/*_Alni_ByOrthoGroup.fasta"): 
		try:
			tmp = str(x)
			tmp= tmp.split("/")
			tmp = tmp[-1]		
			tmp = tmp.split("_")[0]

			mapping_protein_gene = open("/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/transcript_protein_gene/"+ tmp +"_update.txt", "r")
			dict_protein_to_gene = {}
			lines = mapping_protein_gene.readlines()

			for line in lines:
				line = line.split("\t")
				dict_protein_to_gene[line[0]] = line[1]



			mapping_transcript_file = open("tmp/" + tmp + "_mapping_transcript.txt", "r")

			lines = mapping_transcript_file.readlines()
			dict_gene_to_species = {}
			for line in lines:
				line = line.split("\t")
				dict_gene_to_species[line[2]] = line[3]

			
			nt_file = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/seqFastaUsedByEnsembl/" + tmp + ".fasta"
		
			nt_file_2 = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/seqFastaUsedByEnsembl/" + tmp + "_2.fasta"

		
			nt_alignment = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_100_0/" + tmp + "_NT.fasta"
			nt_alignment2 = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_100_0/" + tmp + "_NT_WITH_SPECIES.fasta"
			aa_alignment = "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_100_0/" + tmp + "_AA.fasta"

		
			nt_file_new = ""
			for record in SeqIO.parse(nt_file, "fasta"):
				if record.id in dict_protein_to_gene.keys():
					nt_file_new += ">" + dict_protein_to_gene[record.id] + "\n"
					nt_file_new += str(record.seq) + "\n"

		
			file = open(nt_file_2, "w")
			file.write(nt_file_new)
			file.close()

			compute_alignment_muscle(nt_file_2, x, nt_alignment, aa_alignment)

			with open(nt_alignment) as original, open(nt_alignment2, 'w') as corrected:
				records = SeqIO.parse(nt_alignment, 'fasta')
				for record in records:
					record.id = record.id + "_" + dict_gene_to_species[record.id]
					SeqIO.write(record, corrected, 'fasta')


			command = "./treebest best " + nt_alignment2 + " > " +  "/home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/ensemblGeneTree/sequences_100_0/" + tmp + "_root.nhx -f ressources/ensemblSpecies.tree"
			print (command)
			os.system(command)
		except:
		    pass

recomputeEnsemblTree()
