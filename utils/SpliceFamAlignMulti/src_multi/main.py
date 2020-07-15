#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``main.py`` **module description**:

This module is the main module that 
- parses the input arguments,
- read the input pairwise comparison file
- compute multiple alignment
- computes ortho groups based on multiple alignment
- writes the results in output files

.. moduleauthor:: Aida Ouangraoua

2018

"""


import argparse
import glob
import time
from write import *

from skbio import DistanceMatrix
from skbio.tree import nj

#from compare import *
from compare_progressif import *
from compute_orthology_multi import *


def get_data_from_files(args):
    sourcedata = []
    targetdata = []
    sourcefile = args.sourceFile
    if (sourcefile == None):
        print("Argument -sf <sourcefilename> is required")
    targetfile = args.targetFile
    if (targetfile == None):
        print("Argument -tf <targetfilename> is required")
    source2targetfile = args.source2TargetFile
    if (source2targetfile == None):
        print("Argument -s2tf <source2targetfile> is required")

    if(sourcefile != None and targetfile != None and source2targetfile != None):
        for record in SeqIO.parse(sourcefile, "fasta"):
            sourcedata.append([record.id,str(record.seq),"",[]])
        
        for record in SeqIO.parse(targetfile, "fasta"):
            targetdata.append([record.id,str(record.seq)])
        
        source2target = [line.split("\n")[0].split(" ") for line in open(source2targetfile,"r").readlines()]
        for i in range(len(sourcedata)):
            
            sourcedata[i][2] = source2target[i][1]

        sourceexonfile = args.sourceExonFile
        if(sourceexonfile != None):
            file = open(sourceexonfile, "r")
            lines = file.readlines()
            i = 0
            j = 0
            while j < len(lines):
                line = lines[j]
                if(line.startswith('>')):
                    j += 1
                    line = lines[j]
                    exonlist = []
                    while(j < len(lines) and len(lines[j]) > 2 and not(lines[j].startswith('>'))):
                        line = lines[j]
                        tab = line.split("\n")[0].split(" ")
                        exonlist.append([int(x) for x in tab])
                        j += 1
                    sourcedata[i][3] =  exonlist
                    i += 1
                else:
                    j += 1
    return  sourcedata, targetdata


def parseResultFile(resultfile,sourcedata,idty_threshold):
    file = open(resultfile, "r")
    lines = file.readlines()
    prev_GENE = ""
    liste_GENE    =	[]
    cds2gene = {}
    gene2cds = {}
    cds2geneexon = {}

    for cds in sourcedata:
        cdsid,cdsseq,cdsgeneid,null = cds
        cds2gene[cdsid] = cdsgeneid
        
    comparisonresults			=	[]
    comparisonresults_idty		=	[]
    geneexon = {}
    cdsexon = {}
    liste_CDS_GENE		=	[]
    liste_CDS_GENE_idty		=	[]
    i 					=	0
    while i < len(lines):
        line = lines[i]		
        tab = line.split("\t")
        if len(tab) == 5:
            idCDS  		= 	tab[0]
            idGene 		= 	tab[1]
            i += 1
            CDS_GENE = [0, [], 0, 0, 0]
            CDS_GENE_idty = []
            while i< len(lines) and not(lines[i] in ['\n', '\r\n']):
                line = lines[i]	
                tab = line.split("\t")
                if len(tab) == 9:
                    idCDS			= 	tab[0]
                    idGene			= 	tab[1]
                    length		= 	int(tab[2])
                    beginCDS		= 	int(tab[3])
                    endCDS			= 	int(tab[4])
                    beginGene		=	int(tab[5])
                    endGene			=	int(tab[6])
                    part_CDS_Gene 	= 	[beginCDS, endCDS, beginGene, endGene]
                    part_CDS_Gene_idty = float(tab[7])
                    if(idGene == cds2gene[idCDS]):
                        if(part_CDS_Gene_idty != 1.0):
                            print("WARNING : check exon locations of ",idCDS,"in",idGene, "; exons with less than 100% identity")

                    if(part_CDS_Gene_idty >= idty_threshold and 0 <= beginCDS <= endCDS and 0 <= beginGene <= endGene):
                        CDS_GENE[1].append(part_CDS_Gene)
                        CDS_GENE_idty.append(part_CDS_Gene_idty)
                if(len(tab) == 7):
                    idCDS			= 	tab[0]
                    idGene			= 	tab[1]
                    if(idGene == cds2gene[idCDS]):
                        print("WARNING : check exon locations of ",idCDS,"in",idGene,"; missing exons")

                i += 1

            if(idGene not in geneexon.keys()):
                geneexon[idGene] = []
                gene2cds[idGene]  = []
                
            if(idGene == cds2gene[idCDS]):
                gene2cds[idGene].append(idCDS)
                cdsexon[idCDS] = [x[:2] for x in CDS_GENE[1]]
                geneexon[idGene] += [x[2:] for x in CDS_GENE[1]]
                cds2geneexon[idCDS] = {}
                for x in CDS_GENE[1]:
                    cds2geneexon[idCDS][x[0]] = x[2:]
                
            if(idGene == prev_GENE):
                liste_GENE.append(CDS_GENE)
                liste_GENE_idty.append(CDS_GENE_idty)
            else:
                if(len(liste_GENE) != 0):
                    comparisonresults.append(liste_GENE)
                    comparisonresults_idty.append(liste_GENE_idty)
                liste_GENE    =	[CDS_GENE]
                liste_GENE_idty    =	[CDS_GENE_idty]
            prev_GENE = idGene
        else:
            i += 1

    comparisonresults.append(liste_GENE)
    comparisonresults_idty.append(liste_GENE_idty)
    return comparisonresults,comparisonresults_idty,geneexon,cdsexon,cds2gene,cds2geneexon,gene2cds

def compute_tree(extendedsourcedata,targetdata,comparisonresults,comparisonresults_idty,nbinitialsource):
    rank = {}
    cds2gene = {}
    coverage_matrix = []
    identity_matrix = []
    distance_matrix = []
    list_geneid = []
    for i in range(len(targetdata)):
        coverage_matrix.append([])
        identity_matrix.append([])
        distance_matrix.append([])
        for j in range(len(targetdata)):
            coverage_matrix[i].append([])
            identity_matrix[i].append([])
            distance_matrix[i].append(0)
            
    i=0        
    for gene in targetdata:
        geneid,null = gene
        list_geneid.append(geneid)
        rank[geneid] = i
        i += 1

    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,null,cdsgeneid,null = cds
        cds2gene[cdsid] = cdsgeneid
        j += 1

    i = 0
    for gene in targetdata:
        geneid,null = gene
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            null,blocklist,null,null,null = comparisonresults[i][j]
            coverage = 0
            identity = 0
            k = 0
            for block in blocklist:
                if(block[1] > block[0]):
                    coverage += block[1]-block[0]
                    identity += comparisonresults_idty[i][j][k] * (block[1]-block[0])
                k += 1
            if(coverage != 0):
                identity = 1.0 * identity /coverage
            coverage = 1.0 * coverage /len(cdsseq)
            #print(identity)
            
            coverage_matrix[rank[geneid]][rank[cdsgeneid]].append(coverage)
            coverage_matrix[rank[cdsgeneid]][rank[geneid]].append(coverage)
            identity_matrix[rank[geneid]][rank[cdsgeneid]].append(identity)
            identity_matrix[rank[cdsgeneid]][rank[geneid]].append(identity)
            j += 1
        i += 1
    for i in range(len(targetdata)):
        for j in range(len(targetdata)):
            cover = 0.0
            idty = 0.0
            if(len(coverage_matrix[i][j]) > 0):
               cover = sum(coverage_matrix[i][j])/len(coverage_matrix[i][j])
            if(len(identity_matrix[i][j]) > 0):
               idty = sum(identity_matrix[i][j])/len(identity_matrix[i][j])
            distance_matrix[i][j]= 1.0 - (cover*idty)
            if(i==j and (distance_matrix[i][j] != 0)):
                #print(targetdata[i][0])
                distance_matrix[i][j] = 0
    
    #print(distance_matrix)
    if(len(targetdata) > 3):
        distance_matrix = DistanceMatrix(distance_matrix, list_geneid)
        tree = nj(distance_matrix)
        tree = tree.root_at_midpoint()
        tree = str(tree)
        tree = tree.split("root")[0]+";"
        #print(tree)
    elif(len(targetdata) == 3):
        tree = "(("+list_geneid[0]+","+list_geneid[1]+")"+","+list_geneid[2]+");"
    elif(len(targetdata) == 2):
        tree = "("+list_geneid[0]+","+list_geneid[1]+");"
    else:
        tree = list_geneid[0]+";"
    
    return tree


#####################
### Main ############

def build_arg_parser():
    parser = argparse.ArgumentParser(description="Transcript multiple aligner")
    parser.add_argument('-idty', '--identityThreshold', help="Identity threshold: real between 0.0 and 1.0 (default = 0.3)", type=float, default = 0.3)
    parser.add_argument('-treef', '--treeFile', help="tree file name")
    parser.add_argument('-sf', '--sourceFile', help="Source file name (required)")
    parser.add_argument('-tf', '--targetFile', help="Target file name (required)")
    parser.add_argument('-s2tf', '--source2TargetFile', help="Source to target file name (required)")

    parser.add_argument('-sef', '--sourceExonFile', help="Source exon file name")

    parser.add_argument('-palnf', '--pairwiseAlnFile', help="Pairwise alignment file name (required)")

    parser.add_argument('-op', '--outputPrefix', help="Output prefix (required)")

    parser.add_argument('-ce', '--compareExon', help="The method includes a final step that compares exons of blocks for further merges: Yes or No")

    parser.add_argument('-msa', '--msaMethod', help="Multiple sequence aligner: muscle or mafft")
    
    return parser

def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    treefile = args.treeFile
    idty_threshold = float(args.identityThreshold)
    outputprefix = args.outputPrefix
    pairwisealnfile = args.pairwiseAlnFile

    if (outputprefix == None):
        print("Argument -op <outputPrefix> is required")
        return
            
    if (pairwisealnfile == None):
        print("Argument -palnf <pairwiseAlnFile> is required")
        return

    if(outputprefix != None and pairwisealnfile != None):
        print("Retrieving input data...")
        sourcedata,targetdata = get_data_from_files(args)
        nbinitialsource = len(sourcedata)

        print("Reading pairwise alignment file...")
        comparisonresults,comparisonresults_idty,geneexon,cdsexon,cds2geneid,cds2geneexon,gene2cds = parseResultFile(pairwisealnfile,sourcedata, idty_threshold)
        
        temps=time.time()
        print("Computing multiple alignement...")
        
        tree = ""
        if (treefile != None):
            tree = open(treefile,"r").readline().split("\n")[0]
        else:
            tree = compute_tree(sourcedata,targetdata,comparisonresults,comparisonresults_idty,nbinitialsource)

        #print(tree)
        mblocklist = compute_msa(sourcedata,targetdata,comparisonresults,comparisonresults_idty,geneexon,cdsexon,nbinitialsource,tree,args.compareExon,args.msaMethod,outputprefix)
        print(time.time()-temps)

        temps=time.time()
        print("Writting output files...")
        write_output_files(sourcedata,targetdata,
                           nbinitialsource,mblocklist,outputprefix,args.msaMethod)
        print(time.time()-temps)
        
        temps=time.time()
        print("Computing splicing orthologs...")
        compute_multi_ortholog(cdsexon,outputprefix+"macroalignment.txt", outputprefix+"orthologygrouplist")
        print(time.time()-temps)

if __name__ == '__main__':
    main()
