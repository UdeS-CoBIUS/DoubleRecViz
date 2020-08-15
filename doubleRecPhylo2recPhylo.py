# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : April 2020
@location : University of Sherbrooke

'''
import xml.etree.ElementTree

def dataFromDoubleRecFile(stringio_object):
        recPhylo = []
        spTree = ""
        gnTree = ""
        # file = open(path, "r")
        lines = stringio_object.readlines()
        i = 0
        while(i < len(lines)):
                if("<spTree>" in lines[i]):
                        while("</spTree>" not in lines[i]):
                                spTree += lines[i]
                                i += 1
                        spTree += lines[i]
                        i += 1
                elif("<recGeneTree>" in lines[i]):
                        recPhylo.append([])
                        recGeneTree = ""
                        while("</recGeneTree>" not in lines[i]):
                                recGeneTree += lines[i]
                                i += 1
                        recGeneTree += lines[i]
                        i += 1
                        recPhyloTree = "<recPhylo>\n"+spTree+recGeneTree+"</recPhylo>"
                        recPhylo[-1] = [recPhyloTree,"geneSpecie"]
                        gnTree = recGeneTree2gnTree(recGeneTree)
                elif("<gnTree>" in lines[i]):
                        gnTree = ""
                        while("</gnTree>" not in lines[i]):
                                gnTree += lines[i]
                                i += 1
                        gnTree += lines[i]
                        i += 1
                elif("<recTransTree>" in lines[i]):
                        recPhylo.append([])
                        recTransTree = ""
                        while("</recTransTree>" not in lines[i]):
                                recTransTree += lines[i]
                                i += 1
                        recTransTree += lines[i]
                        i += 1
                        recPhyloTree = "<recPhylo>\n"+gnTree+recTransTree+"</recPhylo>"
                        recPhylo[-1] = [recPhyloTree,"transcriptGene"]
                else:
                        i += 1
        return recPhylo ## array of couples [recPhyloData(string),recType("geneSpecie" or "transcriptGene")]

def recGeneTree2gnTree(recGeneTree):
        gnTree = ""
        recElement = xml.etree.ElementTree.fromstring(recGeneTree)
        root_clade = list(list(recElement)[0])[0]
        gnTree = "<gnTree>\n<phylogeny>\n"+writePrunedTree(root_clade)+"</phylogeny>\n</gnTree>\n"
        return gnTree

def writePrunedTree(clade):
        children = list(clade)
        name = -1
        event = -1
        treeString = ""
        subClades = []
        for i in range(len(children)):
                if children[i].tag == 'name':
                        name = i
                elif children[i].tag == 'eventsRec':
                        event = i
                elif children[i].tag == 'clade':
                        subClades.append(i)
        if(list(children[event])[-1].tag == 'leaf'):
                treeString = "<clade>\n<name>"+children[name].text+"</name>\n</clade>\n"
        else:
                nbSubClades = 0
                for i in subClades:
                        s = writePrunedTree(children[i])
                        if(s != ""):
                               treeString += s
                               nbSubClades += 1
                if(nbSubClades > 1):
                        treeString =  "<clade>\n<name>"+children[name].text+"</name>\n"+treeString+"</clade>\n"
        return treeString
                        
      
