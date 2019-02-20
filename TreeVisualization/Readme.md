# DoubleRecViz
Double reconciliation visualization tool is a web application which allows the double representation of a gene tree embeded inside one specie tree and a transcript tree embeded inside one gene tree.
####Esaie Kuitche, Aïda Ouangraoua, Nadia Tahiri, Jack Parmer and Plotly team
#####Contact Esaie.Kuitche.Kamela@USherbrooke.ca

##Requirements
The program requires the following to be available
- python (2.7 and +) 
- ete3 (https://pypi.python.org/pypi/ete3/)
- cStringIO if you use python2.7 or just IO if you use python3

##Usage
```
python2.7 treeVisualisation.py
```

##Input files
###transcriptGeneReconciliation
Example
```
<recPhylo>
  <gnTree>
    <phylogeny>
      <clade>
        <name>n32</name>
        <clade>
          <name>n30</name>
          <clade>
            <name>n28</name>
            <clade>
              <name>GORILLA</name>
            </clade>
            <clade>
              <name>HUMAN</name>
            </clade>
          </clade>
          <clade>
            <name>MOUSE</name>
          </clade>
        </clade>
        <clade>
          <name>COW</name>
        </clade>
      </clade>
    </phylogeny>
  </gnTree>
  <recTransTree>
    <phylogeny rooted="true">
      <clade>
        <name>n12</name>
        <eventsRec>
          <duplication genesLocation="n32"></duplication>
        </eventsRec>
        <clade>
          <name>n8</name>
          <eventsRec>
            <duplication genesLocation="n30"></duplication>
          </eventsRec>
          <clade>
            <name>r21</name>
            <eventsRec>
              <duplication genesLocation="n28"></duplication>
            </eventsRec>
            <clade>
              <name>LOSS</name>
              <eventsRec>
                <loss genesLocation="GORILLA"></loss>
              </eventsRec>
            </clade>
            <clade>
              <name>r21.c</name>
              <eventsRec>
                <creation genesLocation="HUMAN"></creation>
              </eventsRec>
              <clade>
                <name>gB_human</name>
                <eventsRec>
                  <leaf genesLocation="HUMAN" geneName="gB_human"></leaf>
                </eventsRec>
              </clade>
              <clade>
                <name>gA_human</name>
                <eventsRec>
                  <leaf genesLocation="HUMAN" geneName="gA_human"></leaf>
                </eventsRec>
              </clade>
            </clade>
          </clade>
          <clade>
            <name>n7</name>
            <eventsRec>
              <creation genesLocation="MOUSE"></creation>
            </eventsRec>
            <clade>
              <name>gB_mouse</name>
              <eventsRec>
                <leaf genesLocation="MOUSE"></leaf>
              </eventsRec>
              <clade>
                <name>gA_mouse</name>
                <eventsRec>
                  <leaf genesLocation="MOUSE"></leaf>
                </eventsRec>
              </clade>                      
            </clade>
          </clade>
        </clade>
        <clade>
          <name>cow_all</name>
          <eventsRec>
            <leaf genesLocation="COW"></leaf>
          </eventsRec>                          
        </clade>
      </clade>
    </phylogeny>
  </recTransTree>
</recPhylo>
```
##geneSpecieReconciliation
Example
```
<recPhylo>
  <spTree>
    <phylogeny>
      <clade>
        <name>n32</name>
        <clade>
          <name>n30</name>
          <clade>
            <name>n28</name>
            <clade>
              <name>GORILLA</name>
            </clade>
            <clade>
              <name>HUMAN</name>
            </clade>
          </clade>
          <clade>
            <name>MOUSE</name>
          </clade>
        </clade>
        <clade>
          <name>COW</name>
        </clade>
      </clade>
    </phylogeny>
  </spTree>
  <recGeneTree>
    <phylogeny rooted="true">
      <clade>
        <name>n12</name>
        <eventsRec>
          <speciation speciesLocation="n32"></speciation>
        </eventsRec>
        <clade>
          <name>n8</name>
          <eventsRec>
            <speciation speciesLocation="n30"></speciation>
          </eventsRec>
          <clade>
            <name>r21</name>
            <eventsRec>
              <speciation speciesLocation="n28"></speciation>
            </eventsRec>
            <clade>
              <name>LOSS</name>
              <eventsRec>
                <loss speciesLocation="GORILLA"></loss>
              </eventsRec>
            </clade>
            <clade>
              <name>r21.c</name>
              <eventsRec>
                <duplication speciesLocation="HUMAN"></duplication>
              </eventsRec>
              <clade>
                <name>gB_human</name>
                <eventsRec>
                  <leaf speciesLocation="HUMAN" geneName="gB_human"></leaf>
                </eventsRec>
              </clade>
              <clade>
                <name>gA_human</name>
                <eventsRec>
                  <leaf speciesLocation="HUMAN" geneName="gA_human"></leaf>
                </eventsRec>
              </clade>
            </clade>
          </clade>
          <clade>
            <name>n7</name>
            <eventsRec>
              <duplication speciesLocation="MOUSE"></duplication>
            </eventsRec>
            <clade>
              <name>gB_mouse</name>
              <eventsRec>
                <leaf speciesLocation="MOUSE"></leaf>
              </eventsRec>
              <clade>
                <name>gA_mouse</name>
                <eventsRec>
                  <leaf speciesLocation="MOUSE"></leaf>
                </eventsRec>
              </clade>                      
            </clade>
          </clade>
        </clade>
        <clade>
          <name>cow_all</name>
          <eventsRec>
            <leaf speciesLocation="COW"></leaf>
          </eventsRec>                          
        </clade>
      </clade>
    </phylogeny>
  </recGeneTree>
</recPhylo>
```

outfile name
```
the reconciliated trees are draw in your browser
```
##Running DoubleRecViz on an example
```
python2.7 treeVisualisation.py
```
