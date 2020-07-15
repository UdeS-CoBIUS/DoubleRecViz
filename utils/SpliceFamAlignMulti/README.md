# SpliceFamAlignMulti
Program Version 2019

Performs a multiple spliced alignment from the set of all pairs of source CDS and target gene comparison in a progressive way.
----------------------------------------------------------------

Authors: Abigail Djossou, Safa Jammali, and Aida Ouangraoua

Université de Sherbrooke, Canada
CoBIUS Lab:  https://cobius.usherbrooke.ca/

For questions email us at abigail.djossou@usherbrooke.ca

### Requirements:

-Bio
-Mafft+ standalone
-Muscle+ standalone
-Blast+ standalone
-collections
-newick
-numpy
-contextlib
-functools
-multiprocessing
-argparse


### Usage
```
pairwise: python2.7 main.py [-h] [-c CHOICESTRUCTURE] [-s STEP] [-ce COMPAREEXON]
               [-sf SOURCEFILE] [-tf TARGETFILE] [-s2tf SOURCE2TARGETFILE]
               [-sef SOURCEEXONFILE] [-op OUTPUTPREFIX] [-of OUTPUTFORMAT]

```
*  -h, --help            show this help message and exit
*  -c CHOICESTRUCTURE, --choiceStructure CHOICESTRUCTURE
                        Method used to infer splicing structure when a
                        splicing structure file is not given: blast or splign
*  -s STEP, --step STEP  The method goes until Step 1, 2, or 3: 1, 2 or 3
                        (required)
*  -ce COMPAREEXON, --compareExon COMPAREEXON
                        The method includes in Step2 a comparison of exons:
                        Yes or No (required)
*  -sf SOURCEFILE, --sourceFile SOURCEFILE
                        Source (CDS) file name (required)
*  -tf TARGETFILE, --targetFile TARGETFILE
                        Target (gene) file name (required)
*  -s2tf SOURCE2TARGETFILE, --source2TargetFile SOURCE2TARGETFILE
                        Association between source and target file name
                        (required)
*  -sef SOURCEEXONFILE, --sourceExonFile SOURCEEXONFILE
                        Source Exon (splicing structure) file name
*  -op OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Output prefix (required)
*  -of OUTPUTFORMAT, --outputFormat OUTPUTFORMAT
                        Output format : list or aln (required)   

```
multiple: python3 main.py [-h] [-idty IDENTITYTHRESHOLD] [-treef TREEFILE]
               [-sf SOURCEFILE] [-tf TARGETFILE] [-s2tf SOURCE2TARGETFILE]
               [-sef SOURCEEXONFILE] [-palnf PAIRWISEALNFILE]
               [-op OUTPUTPREFIX] [-ce COMPAREEXON] [-msa MSAMETHOD]

```
*  -h, --help            show this help message and exit
*  -idty IDENTITYTHRESHOLD, --identityThreshold IDENTITYTHRESHOLD
                        Identity threshold: real between 0.0 and 1.0 (default
                        = 0.3)
*  -treef TREEFILE, --treeFile TREEFILE
                        tree file name
*  -sf SOURCEFILE, --sourceFile SOURCEFILE
                        Source file name (required)
*  -tf TARGETFILE, --targetFile TARGETFILE
                        Target file name (required)
*  -s2tf SOURCE2TARGETFILE, --source2TargetFile SOURCE2TARGETFILE
                        Source to target file name (required)
*  -sef SOURCEEXONFILE, --sourceExonFile SOURCEEXONFILE
                        Source exon file name
*  -palnf PAIRWISEALNFILE, --pairwiseAlnFile PAIRWISEALNFILE
                        Pairwise alignment file name (required)
*  -op OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Output prefix (required)
*  -ce COMPAREEXON, --compareExon COMPAREEXON
                        The method includes a final step that compares exons
                        of blocks for further merges: Yes or No
*  -msa MSAMETHOD, --msaMethod MSAMETHOD
                        Multiple sequence aligner: muscle or mafft
### Running SpliceFamAlignMulti: examples of command line

#### First, compute all pairwise alignments (Python2.7):
```
python2.7 src/main.py -s 2 -ce Yes -sf examples/input/FAM86SG/FAM86SG_initialsource.fasta -tf examples/input/FAM86SG/FAM86SG_target.fasta -s2tf examples/input/FAM86SG/FAM86SG_initialsource2target.txt -sef examples/input/FAM86SG/FAM86SG_initialsourceexonlist.txt -op examples/output/FAM86SG_ -of list
```
#### Then, compute multiple spliced alignment  (Python3):
```
python3 src_multi/main.py -ce Yes -sf examples/input/FAM86SG/FAM86SG_initialsource.fasta -tf examples/input/FAM86SG/FAM86SG_target.fasta -s2tf examples/input/FAM86SG/FAM86SG_initialsource2target.txt -sef examples/input/FAM86SG/FAM86SG_initialsourceexonlist.txt -palnf examples/output/FAM86SG_result.txt -op examples/output/FAM86SG_ 
```
