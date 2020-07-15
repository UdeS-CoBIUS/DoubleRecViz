/!\ Les chemins sont écrits en durs, ils sont à changer dans le code pour qu'ils marchent sur votre machine 
Liste des fichiers dans lesquels il faut changer les chemins:
spliceGraphMaker.py
alternativeSplicing.py
labelExons.py
matrixScore.py
clustering.py
retrieveSeqClusters.py
retrieveSeqClustersDir.py
microWithoutGenes.py
fusionMatrix.py 
diffmatrix.py
/!\ 

##Labelisation des exons :
- AVEC SPLICEGRAPHE 
Génère les splicegraphes pour le répertoire macroalignement, résultats dans le répertoire spliceGraph : 
<<<<<<< HEAD
python spliceGraphMaker.py

Labelisation des exons, résultats dans le répertoire labelExonsMaker : 

 python alternativeSplicing.py -i spliceGraph/FAM86_spliceGraphe.txt -ma macroalignment/FAM86_macroalignment.txt -mi microalignment/FAM86_microalignment.fasta
=======

- AVEC MICROALIGNEMENT
Labelisation des exons, résultats dans le répertoire labelExons : 

python labelExons.py macroalignment/FAM86_macroalignment.txt microalignment/FAM86_microalignment.fasta


Matrice de similarités :
Résultats dans similarityScores : 

Commande avec le microalignment
python -W ignore matrixScore.py labelExons/FAM86_labeledExons.out initialSource/FAM86_initialsource2target.txt


Commande avec le spliceGraph
python -W ignore matrixScore.py labelExonsMaker/FAM86_spliceGraphe_labeledExons.out initialSource/FAM86_initialsource2target.txt

Pour lancer le script sur le répertoire labelExons ou labelExonsMaker : 

python scoreDirectory.py
Ne pas oublier de donner les droits d'éxécution au script : chmod +x matrixScore.py


Clustering ascendant hiérarchique : 
Résultats dans le répertoire clusters: 

python -W ignore clustering.py similarityScores/FAM86_score.csv 0.5
Le deuxième argument est le minimum de similarités qu'on veut entre deux transcripts

Pour le lancer sur un répertoire : python scoreDirectory.py
Ne pas oublier de donner les droits d'éxécution au script : chmod +x clustering.py

Récupération des séquences du microalignements pour chaque sous groupes : 
Résultats dans sequencesClusters, chaque famille a son propre répertoire avec ses sous groupes : 

python retrieveSeqClusters.py clusters/FAM86_clusters.out microaligment/FAM86_microalignment.fasta

Pour lancer le script sur un répertoire : python retrieveSeqClustersDir.py
Ne pas oublier de donner les droits d'éxécution au script :  chmod +x retrieveSeqClusters.py

Conversion des fichiers des séquences des sous groupes (.fasta) en fichier phylip : 
Résultats dans le répertoire phylipFormat d'une famille donnée : 
python fasta_to_phylip.py seqs.fasta seqs.phylip


Pour exécuter le script sur tout les sous groupes de toutes les familles : python fasta_to_phylip_dir.py sequenceClusters_directory_path

Pour obtenir les arbres des sous groupes, il faut utiliser PhyML : 
Pour l'installer dans le bash du home :
$ nano ~/.bashrc
export PATH="/home/local/USHERBROOKE/degm2303/Documents/SpliceGraph/PhyML-3.1/:$PATH
Lancer le script python phyml_tree.py /home/local/USHERBROOKE/degm2303/Documents/SpliceGraph/clusters/sequencesClusters qui créera l'arbre pour chaque sous groupe de chaque famille du répertoire donné

Enraciner les arbres grâce à Njplot : ./njplot


Fusion des matrices de structures et des matrices de séquences : 
Obtenir la matrice de séquences, dans le dossier scripts/FSePSA (il est au meme niveau que le dossier SpliceGraph)
Lancer le script:

cd ../scripts/FSePSA/examples 

python microWithoutGenes.py FAM86_microalignment.fasta

pour enlever les genes du microalignements

cd ../src/

python compute_score_main.py 
(changer les paramètres --dataalignment et --outfile pour correspondre au nom de la famille)


cp /home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/scripts/FsePSA/examples/FAM86.matrix.txt /home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph/similaritySeqScores/

cd /home/local/USHERBROOKE/kuie2201/Bureau/MarieDegen/SpliceGraph

python -W ignore fusionMatrix.py similaritySeqScores/FAM86.matrix.txt similarityScores/FAM86_score.csv

Résultats dans le répertoire scores/ 


Différences des deux matrices obtenues a partir des deux methodes:
python -W ignore diffmatrix.py similarityScores/FAM86_score.csv similarityScores/FAM86_score1.csv



