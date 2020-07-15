#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``main.py`` **module description**:

This module is the main module that 
- parses the input arguments,
- launches the comparison for all pairs of source CDS and target gene,
- computes splicing ortholog groups
- writes the results in output files

Command line: 
If a splicing structure file is not given in input:
python 'pathto main.py'  -sf 'path to cds fasta file'  -tf 'path to gene fasta file' -s2tf 'path to association cds-gene file'  -op 'output path' -of 'list or aln: display alignments as lists of blocks or detailed alignments'  -f 'Yes or No: if the method goes until the SFA_G step or not' -c 'blast or splign: method used to infer splicing structure when a splicing structure file is not given'

If a splicing structure file is given in input:
python 'pathto main.py'  -sf 'path to cds fasta file'  -tf 'path to gene fasta file' -s2tf 'path to association cds-gene file' -sef 'splicing structure file' -op 'output path' -of 'list or aln: display alignments as lists of blocks or detailed alignments -f 'Yes or No: if the method goes until the SFA_G step or not'
 

moduleauthor:: Safa Jammali, Jean-David Aguilar and Aïda Ouangraoua 
Université de Sherbrooke, Canada
for questions email us at Safa.Jammali@USherbrooke.ca

2017-2018

"""

# import external bilio

import os
import argparse
import glob
import time

# import project-specific file
from read import *
from compute_alignments import *
from compute_orthologs import *
from write import *



#####################
### Main ############


def build_arg_parser():
    """
    This function parses all parameters given by the user in th command line

    Parameters
    ----------


    Returns
    -------
    parser
        list of arguments
    """

    parser = argparse.ArgumentParser(description="SpliceFamAlign")
    parser.add_argument('-c', '--choiceStructure', help="Method used to infer splicing structure when a splicing structure file is not given: blast or splign", default = "splign")
    parser.add_argument('-pm', '--pairwiseMethod', help="Method used to compute pairwise spliced alignment: sfa or splign", default = "sfa")
    parser.add_argument('-s', '--step', help="The method goes until Step 1, 2, or 3: 1, 2 or 3 (required)")
    parser.add_argument('-ce', '--compareExon', help="The method includes in Step2 a comparison of exons: Yes or No (required)")
    
    #required for it = file"
    parser.add_argument('-sf', '--sourceFile', help="Source (CDS) file name (required)")
    parser.add_argument('-tf', '--targetFile', help="Target (gene) file name (required)")
    parser.add_argument('-s2tf', '--source2TargetFile', help="Association between source and  target file name (required)")
    parser.add_argument('-sef', '--sourceExonFile', help="Source Exon (splicing structure) file name")
            
    parser.add_argument('-op', '--outputPrefix', help="Output prefix (required)")

    parser.add_argument('-of', '--outputFormat', help="Output format : list or aln (required)")

    return parser

def main():
    """
    This is the main function, it takes all input arguments and launches
        all module of the project
        
    Parameters
    ----------


    Returns
    -------
    
    """
    files1=glob.glob('./sequences/genes/*')
    for file1 in files1:
        os.remove(file1)
    files2=glob.glob('./sequences/cds/*')
    for file2 in files2:
        os.remove(file2)
    files3 = glob.glob('./results/splign_results/*')
    for file3 in files3:
        os.remove(file3)
    files4 = glob.glob('./results/ident_results/*')
    for file4 in files4:
        os.remove(file4)
    files5 = glob.glob('./results/blast_Results/*')
    for file5 in files5:
        os.remove(file5)
    # parses the input arguments
    parser = build_arg_parser()
    args = parser.parse_args()
    outputPrefix = args.outputPrefix
    outputFormat = args.outputFormat
    choice = args.choiceStructure
    pairwiseMethod = args.pairwiseMethod
    step = args.step
    compareExon = args.compareExon
    
    if (outputPrefix == None):
        print "Argument -op <outputprefix> is required"

    if (outputFormat == None):
        print "Argument -of <outputformat> is required : list or aln"

    if(outputPrefix != None and outputFormat != None):
        print "Retrieving input data..."
        sourceData,targetData = get_data_from_files(args)
        nbinitialSource = len(sourceData)

        # Compute all  pairwise alignments for source CDS and target gene
        temps=time.time()
        print "Comparing sequences..."
        comparisonResults = spliceAlignment(sourceData, targetData, step, compareExon, choice, pairwiseMethod, outputPrefix)
        print(time.time()-temps)

        # Compute splicing orthologs
        temps=time.time()
        print "Computing splicing ortholog groups..."
        orthologyGroups = computeOrthology(sourceData,targetData,comparisonResults)
        print(time.time()-temps)

        temps=time.time()
        print "Completing orthology groups with putative non-annotated CDS..."
        extendedSourceData,extendedOrthologyGroups = completeOrthology(sourceData,targetData,comparisonResults,orthologyGroups)
        print(time.time()-temps)
        
        # write the results in output files
        temps=time.time()
        print "Writting output files..."        
        writeOutfile(outputPrefix,outputFormat,extendedSourceData,targetData,
                     extendedOrthologyGroups,
                     comparisonResults,nbinitialSource)
        print(time.time()-temps)

        
if __name__ == '__main__':
    main()

