#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``getdata.py`` **module description**:

This module computes the input data structures using the input files

.. moduleauthor::  Safa Jammali and  Aida Ouangraoua


2017-2018

"""

from Bio import SeqIO


#########################
## DATA ACQUISITION #####
#########################

##### get data from files ##########

def get_data_from_files(args):
    sourcedata = []
    targetdata = []
    sourcefile = args.sourceFile
    if (sourcefile == None):
        print "Argument -sf <sourcefilename> is required"
    targetfile = args.targetFile
    if (targetfile == None):
        print "Argument -tf <targetfilename> is required"
    source2targetfile = args.source2TargetFile
    if (source2targetfile == None):
        print "Argument -s2tf <source2targetfile> is required"

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
                    while( j < len(lines) and len(lines[j]) > 2 and not(lines[j].startswith('>'))):
                        line = lines[j]
                        tab = line.split("\n")[0].split(" ")
                        exonlist.append([int(x) for x in tab])
                        j += 1
                    sourcedata[i][3] =  exonlist
                    i+=1
                else:
                    j += 1
    return  sourcedata, targetdata



