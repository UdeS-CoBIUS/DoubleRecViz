#!/bin/bash
import os

# Pairwise ENSGT00390000009061
cmd1 = "python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00390000009061/ENSGT00390000009061_initialsource.fasta -tf examples/input/ENSGT00390000009061/ENSGT00390000009061_target.fasta -s2tf examples/input/ENSGT00390000009061/ENSGT00390000009061_initialsource2target.txt -sef examples/input/ENSGT00390000009061/ENSGT00390000009061_initialsourceexonlist.txt -op examples/output/ENSGT00390000009061_ -of list"

# Multi ENSGT00390000009061
cmd2 = "python3 src_multi/main.py -sf examples/input/ENSGT00390000009061/ENSGT00390000009061_initialsource.fasta -tf examples/input/ENSGT00390000009061/ENSGT00390000009061_target.fasta -s2tf examples/input/ENSGT00390000009061/ENSGT00390000009061_initialsource2target.txt -sef examples/input/ENSGT00390000009061/ENSGT00390000009061_initialsourceexonlist.txt -palnf examples/output/ENSGT00390000009061_result.txt -op examples/output/ENSGT00390000009061_ -ce Yes"

os.system(cmd1)
os.system(cmd2)
