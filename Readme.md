# DoubleRecViz
DoubleRecViz is a tool for visualizing double reconciliations between phylogenetic trees at three levels: transcript, gene and species. It displays reconciled gene trees embedded in species trees and reconciled transcript trees embedded in gene trees, with annotations on the nodes of the trees corresponding to the underlying co-evolutionary events.

#### Esaie Kuitche, Yanchun Qi, Nadia Tahiri, Jack Parmer, Aïda Ouangraoua

##### Contact Esaie.Kuitche.Kamela@USherbrooke.ca

## Requirements
The program requires the following to be available
- Plotly and Dash libraries
- Python (2.7 and +) 
- ETE3 
- cStringIO if you use python2.7 or just IO if you use python3

A procedure for the installation of requirements can be found in file "Installation.txt".

## Usage
```
python treeVisualisation.py
```

## Input files

DoubleRecViz makes use of the doubleRecPhyloXML and recTransTreeXML grammars which are extensions of the recPhyloXML and recGeneTreeXML grammars inherited from phyloXML and designed to describe reconciled gene-species tree reconciliations (see https://github.com/WandrilleD/recPhyloXML/tree/master/xsd for a detailed description of the recPhyloXML and recGeneTreeXML grammars). A detailed description of the doubleRecPhyloXML and recTransTreeXML grammars can be found in "manual.html" and in the "XSD" directory. 

DoubleRecViz takes as input one the following types of object:

 - a doubleRecPhylo object containing 1 species tree (spTree object) followed by 1 to n sets such that each set is composed of 1 reconciled gene tree (recGeneTree object) followed by 0 to n reconciled transcript trees (recTranscriptTree object).

 - a recPhylo object containing 1 species tree (spTree object) followed by 1 to n reconciled gene trees (recGeneTree object)

 - a recPhylo object containing 1 gene tree (gnTree object) followed by 1 to n reconciled transcript trees (recTransTree object) 

Examples input file data (doubleRecPhylo objects for a transcript-gene-species tree reconciliation, and recPhylo objects for a transcript-gene tree reconciliations and gene-species tree reconciliations) can be found on the "testData" directory.

## Output

```
The reconciled trees are displayed in the browser.
```

## Generation of test data from sequence data


## Conversion between the NHX and doubleRecPhyloXML formats

Utility scripts for conversion between the NHX and doubleRecPhyloXML formats are also provided in the "convertFormat" directory. Use -h option when running the scripts to display the description.

Examples of commands to convert from  NHX to doubleRecPhyloXML:

```
python convertFormat/NHXtoDoubleRecPhyloXML.py -t double -s testData/realData/ENSGT00390000015814/Trees/ENSGT00390000015814_pruned_speciestree.nw -rgt testData/realData/ENSGT00390000015814_double.nhx 
```

```
python convertFormat/NHXtoDoubleRecPhyloXML.py -t genespecies -s testData/realData/ENSGT00390000015814/Trees/ENSGT00390000015814_pruned_speciestree.nw -rg testData/realData/ENSGT00390000015814_genespecies.nhx 
```

```
python convertFormat/NHXtoDoubleRecPhyloXML.py -t transcriptgene -g testData/realData/ENSGT00390000015814/Trees/ENSGT00390000015814_genetree_init.nw -rt testData/realData/ENSGT00390000015814_transcriptgene.nhx
```
