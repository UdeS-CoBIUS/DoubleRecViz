#!/bin/bash

# Pairwise FAM86SG
python src/main.py -s 2 -ce Yes -sf examples/input/FAM86SG/FAM86SG_initialsource.fasta -tf examples/input/FAM86SG/FAM86SG_target.fasta -s2tf examples/input/FAM86SG/FAM86SG_initialsource2target.txt -sef examples/input/FAM86SG/FAM86SG_initialsourceexonlist.txt -op examples/output/FAM86SG_ -of list

# Multi FAM86SG
python3 src_multi/main.py -sf examples/input/FAM86SG/FAM86SG_initialsource.fasta -tf examples/input/FAM86SG/FAM86SG_target.fasta -s2tf examples/input/FAM86SG/FAM86SG_initialsource2target.txt -sef examples/input/FAM86SG/FAM86SG_initialsourceexonlist.txt -palnf examples/output/FAM86SG_result.txt -op examples/output/FAM86SG_ -ce Yes

# Pairwise FAM86
python src/main.py -s 2 -ce Yes -sf examples/input/FAM86/FAM86_initialsource.fasta -tf examples/input/FAM86/FAM86_target.fasta -s2tf examples/input/FAM86/FAM86_initialsource2target.txt -sef examples/input/FAM86/FAM86_initialsourceexonlist.txt -op examples/output/FAM86_ -of list

# Multi FAM86
python3 src_multi/main.py -sf examples/input/FAM86/FAM86_initialsource.fasta -tf examples/input/FAM86/FAM86_target.fasta -s2tf examples/input/FAM86/FAM86_initialsource2target.txt -sef examples/input/FAM86/FAM86_initialsourceexonlist.txt -palnf examples/output/FAM86_result.txt -op examples/output/FAM86_ -ce Yes

# Pairwise MAG
python src/main.py -s 2 -ce Yes -sf examples/input/MAG/MAG_initialsource.fasta -tf examples/input/MAG/MAG_target.fasta -s2tf examples/input/MAG/MAG_initialsource2target.txt -sef examples/input/MAG/MAG_initialsourceexonlist.txt -op examples/output/MAG_ -of list

# Multi MAG
python3 src_multi/main.py -sf examples/input/MAG/MAG_initialsource.fasta -tf examples/input/MAG/MAG_target.fasta -s2tf examples/input/MAG/MAG_initialsource2target.txt -sef examples/input/MAG/MAG_initialsourceexonlist.txt -palnf examples/output/MAG_result.txt -op examples/output/MAG_ -ce Yes

# Pairwise ENSGT00940000158754
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource.fasta -tf examples/input/ENSGT00940000158754/ENSGT00940000158754_target.fasta -s2tf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource2target.txt -sef examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsourceexonlist.txt -op examples/output/ENSGT00940000158754_ -of list

# Multi ENSGT00940000158754
python3 src_multi/main.py -sf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource.fasta -tf examples/input/ENSGT00940000158754/ENSGT00940000158754_target.fasta -s2tf examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsource2target.txt -sef examples/input/ENSGT00940000158754/ENSGT00940000158754_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000158754_result.txt -op examples/output/ENSGT00940000158754_ -ce Yes

# Pairwise ENSGT00940000162816
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource.fasta -tf examples/input/ENSGT00940000162816/ENSGT00940000162816_target.fasta -s2tf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource2target.txt -sef examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsourceexonlist.txt -op examples/output/ENSGT00940000162816_ -of list

# Multi ENSGT00940000162816
python3 src_multi/main.py -sf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource.fasta -tf examples/input/ENSGT00940000162816/ENSGT00940000162816_target.fasta -s2tf examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsource2target.txt -sef examples/input/ENSGT00940000162816/ENSGT00940000162816_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000162816_result.txt -op examples/output/ENSGT00940000162816_  -ce Yes

# Pairwise ENSGT00390000005420
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00390000005420/ENSGT00390000005420_initialsource.fasta -tf examples/input/ENSGT00390000005420/ENSGT00390000005420_target.fasta -s2tf examples/input/ENSGT00390000005420/ENSGT00390000005420_initialsource2target.txt -sef examples/input/ENSGT00390000005420/ENSGT00390000005420_initialsourceexonlist.txt -op examples/output/ENSGT00390000005420_ -of list

# Multi ENSGT00390000005420
python3 src_multi/main.py -sf examples/input/ENSGT00390000005420/ENSGT00390000005420_initialsource.fasta -tf examples/input/ENSGT00390000005420/ENSGT00390000005420_target.fasta -s2tf examples/input/ENSGT00390000005420/ENSGT00390000005420_initialsource2target.txt -sef examples/input/ENSGT00390000005420/ENSGT00390000005420_initialsourceexonlist.txt -palnf examples/output/ENSGT00390000005420_result.txt -op examples/output/ENSGT00390000005420_  -ce Yes

# Pairwise ENSGT00910000148697
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00910000148697/ENSGT00910000148697_initialsource.fasta -tf examples/input/ENSGT00910000148697/ENSGT00910000148697_target.fasta -s2tf examples/input/ENSGT00910000148697/ENSGT00910000148697_initialsource2target.txt -sef examples/input/ENSGT00910000148697/ENSGT00910000148697_initialsourceexonlist.txt -op examples/output/ENSGT00910000148697_ -of list

# Multi ENSGT00910000148697
python3 src_multi/main.py -sf examples/input/ENSGT00910000148697/ENSGT00910000148697_initialsource.fasta -tf examples/input/ENSGT00910000148697/ENSGT00910000148697_target.fasta -s2tf examples/input/ENSGT00910000148697/ENSGT00910000148697_initialsource2target.txt -sef examples/input/ENSGT00910000148697/ENSGT00910000148697_initialsourceexonlist.txt -palnf examples/output/ENSGT00910000148697_result.txt -op examples/output/ENSGT00910000148697_  -ce Yes

# Pairwise ENSGT00940000166277
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000166277/ENSGT00940000166277_initialsource.fasta -tf examples/input/ENSGT00940000166277/ENSGT00940000166277_target.fasta -s2tf examples/input/ENSGT00940000166277/ENSGT00940000166277_initialsource2target.txt -sef examples/input/ENSGT00940000166277/ENSGT00940000166277_initialsourceexonlist.txt -op examples/output/ENSGT00940000166277_ -of list

# Multi ENSGT00940000166277
python3 src_multi/main.py -sf examples/input/ENSGT00940000166277/ENSGT00940000166277_initialsource.fasta -tf examples/input/ENSGT00940000166277/ENSGT00940000166277_target.fasta -s2tf examples/input/ENSGT00940000166277/ENSGT00940000166277_initialsource2target.txt -sef examples/input/ENSGT00940000166277/ENSGT00940000166277_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000166277_result.txt -op examples/output/ENSGT00940000166277_  -ce Yes

# Pairwise ENSGT00940000167924
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000167924/ENSGT00940000167924_initialsource.fasta -tf examples/input/ENSGT00940000167924/ENSGT00940000167924_target.fasta -s2tf examples/input/ENSGT00940000167924/ENSGT00940000167924_initialsource2target.txt -sef examples/input/ENSGT00940000167924/ENSGT00940000167924_initialsourceexonlist.txt -op examples/output/ENSGT00940000167924_ -of list

# Multi ENSGT00940000167924
python3 src_multi/main.py -sf examples/input/ENSGT00940000167924/ENSGT00940000167924_initialsource.fasta -tf examples/input/ENSGT00940000167924/ENSGT00940000167924_target.fasta -s2tf examples/input/ENSGT00940000167924/ENSGT00940000167924_initialsource2target.txt -sef examples/input/ENSGT00940000167924/ENSGT00940000167924_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000167924_result.txt -op examples/output/ENSGT00940000167924_  -ce Yes

# Pairwise ENSGT00940000178467
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00940000178467/ENSGT00940000178467_initialsource.fasta -tf examples/input/ENSGT00940000178467/ENSGT00940000178467_target.fasta -s2tf examples/input/ENSGT00940000178467/ENSGT00940000178467_initialsource2target.txt -sef examples/input/ENSGT00940000178467/ENSGT00940000178467_initialsourceexonlist.txt -op examples/output/ENSGT00940000178467_ -of list

# Multi ENSGT00940000178467
python3 src_multi/main.py -sf examples/input/ENSGT00940000178467/ENSGT00940000178467_initialsource.fasta -tf examples/input/ENSGT00940000178467/ENSGT00940000178467_target.fasta -s2tf examples/input/ENSGT00940000178467/ENSGT00940000178467_initialsource2target.txt -sef examples/input/ENSGT00940000178467/ENSGT00940000178467_initialsourceexonlist.txt -palnf examples/output/ENSGT00940000178467_result.txt -op examples/output/ENSGT00940000178467_  -ce Yes

# Pairwise ENSGT00950000182910
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00950000182910/ENSGT00950000182910_initialsource.fasta -tf examples/input/ENSGT00950000182910/ENSGT00950000182910_target.fasta -s2tf examples/input/ENSGT00950000182910/ENSGT00950000182910_initialsource2target.txt -sef examples/input/ENSGT00950000182910/ENSGT00950000182910_initialsourceexonlist.txt -op examples/output/ENSGT00950000182910_ -of list

# Multi ENSGT00950000182910
python3 src_multi/main.py -sf examples/input/ENSGT00950000182910/ENSGT00950000182910_initialsource.fasta -tf examples/input/ENSGT00950000182910/ENSGT00950000182910_target.fasta -s2tf examples/input/ENSGT00950000182910/ENSGT00950000182910_initialsource2target.txt -sef examples/input/ENSGT00950000182910/ENSGT00950000182910_initialsourceexonlist.txt -palnf examples/output/ENSGT00950000182910_result.txt -op examples/output/ENSGT00950000182910_  -ce Yes

# Pairwise ENSGT00960000189321
python src/main.py -s 2 -ce Yes -sf examples/input/ENSGT00960000189321/ENSGT00960000189321_initialsource.fasta -tf examples/input/ENSGT00960000189321/ENSGT00960000189321_target.fasta -s2tf examples/input/ENSGT00960000189321/ENSGT00960000189321_initialsource2target.txt -sef examples/input/ENSGT00960000189321/ENSGT00960000189321_initialsourceexonlist.txt -op examples/output/ENSGT00960000189321_ -of list

# Multi ENSGT00960000189321
python3 src_multi/main.py -sf examples/input/ENSGT00960000189321/ENSGT00960000189321_initialsource.fasta -tf examples/input/ENSGT00960000189321/ENSGT00960000189321_target.fasta -s2tf examples/input/ENSGT00960000189321/ENSGT00960000189321_initialsource2target.txt -sef examples/input/ENSGT00960000189321/ENSGT00960000189321_initialsourceexonlist.txt -palnf examples/output/ENSGT00960000189321_result.txt -op examples/output/ENSGT00960000189321_  -ce Yes

