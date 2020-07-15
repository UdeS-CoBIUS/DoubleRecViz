#!/bin/bash

cd SpliceFamAlignMulti
# Pairwise ENSGT00390000015814
python2.7 src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00390000015814/ENSGT00390000015814_initialsource.fasta -tf examples/input/ENSGT00390000015814/ENSGT00390000015814_target.fasta -s2tf examples/input/ENSGT00390000015814/ENSGT00390000015814_initialsource2target.txt -sef examples/input/ENSGT00390000015814/ENSGT00390000015814_initialsourceexonlist.txt -op examples/output/ENSGT00390000015814_ -of list

# Multi ENSGT00390000015814
python2.7 src_multi/main.py -sf examples/input/ENSGT00390000015814/ENSGT00390000015814_initialsource.fasta -tf examples/input/ENSGT00390000015814/ENSGT00390000015814_target.fasta -s2tf examples/input/ENSGT00390000015814/ENSGT00390000015814_initialsource2target.txt -sef examples/input/ENSGT00390000015814/ENSGT00390000015814_initialsourceexonlist.txt -palnf examples/output/ENSGT00390000015814_result.txt -op examples/output/ENSGT00390000015814_ -ce Yes

# moves files to SHT program
cd ..
python2.7 moveFiles.py -g ENSGT00390000015814

#run SHT program
cd SHT
python2.7 spliceGraphMaker.py -i macroalignment/ENSGT00390000015814_macroalignment.txt

#Extract SHT gene tree, and set of transcript trees
cd ..
python2.7  extractTree.py -g ENSGT00390000015814

#Run SuperProteineTree to combine all transcript trees
cd SuperProteinTree
python2.7 src/main.py -i ENSGT00390000015814

#Compute Double Reconciliation between transcript, gene and species tree.
cd ..
python2.7 computeDoubleReconciliation.py -id ENSGT00390000015814

#Move result files in its directory
python2.7 writeResult.py -id ENSGT00390000015814
