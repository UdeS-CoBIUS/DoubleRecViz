# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Esaie Kuitche
@date : April 2020

@location : University of Sherbrooke

'''
import xml.etree.ElementTree
from Bio import Phylo
from ete3 import Tree
from math import *

def make_dict_from_tree(element_tree, recType):
    """Traverse the given XML element tree to convert it into a dictionary.
    https://ericscrivner.me/2015/07/python-tip-convert-xml-tree-to-a-dictionary/ 
    :param element_tree: An XML element tree
    :type element_tree: xml.etree.ElementTree
    :rtype: dict
    """
    def internal_iter(tree, accum, recType):
        """Recursively iterate through the elements of the tree accumulating
        a dictionary result.
 
        :param tree: The XML element tree
        :type tree: xml.etree.ElementTree
        :param accum: Dictionary into which data is accumulated
        :type accum: dict
        :rtype: dict
        """
        if tree is None:
            return accum
 
        if tree.getchildren():
            accum[tree.tag] = {}
            for each in tree.getchildren():
                result = internal_iter(each, {}, recType)
                if each.tag in accum[tree.tag]:
                    if not isinstance(accum[tree.tag][each.tag], list):
                        accum[tree.tag][each.tag] = [
                            accum[tree.tag][each.tag]
                        ]
                    accum[tree.tag][each.tag].append(result[each.tag])
                else:
                    accum[tree.tag].update(result)
        else:
            if tree.tag in ["speciation", "duplication", "branchingOut", "loss", "leaf", "creation"]:
                if recType == "geneSpecie":
                    accum[tree.tag] = {"speciesLocation":tree.attrib['speciesLocation']}
                elif recType == "transcriptGene":
                    accum[tree.tag] = {"genesLocation":tree.attrib['genesLocation']}                
            elif tree.tag in ["bifurcationOut"]:
                accum[tree.tag] = {"bifurcationOut":tree.text}
            elif tree.tag in ["transferBack"]:
                accum[tree.tag] =  {"transferBack":tree.attrib['destinationSpecies']}
            else:
                accum[tree.tag] = tree.text
 
        return accum
 
    return internal_iter(element_tree, {}, recType)
 
 
def buildTreeFromXml(sptre, tree):
    if type(sptre) ==  dict:
        if  len(sptre.keys()) == 1:
            k = list(sptre.keys())[0]
            v = sptre[k]
            if len(v) == 2:
                tree.add_child(name=v['name'])
                clade = tree.search_nodes(name=v["name"])[0]
                return buildTreeFromXml(v["clade"], clade)
            else:
                tree.add_child(name=sptre['name'])
                return tree
                
        elif  len(sptre.keys()) == 2:
                tree.add_child(name=sptre['name'])
                clade = tree.search_nodes(name=sptre["name"])[0]
                return buildTreeFromXml(sptre["clade"], clade)        
    elif type(sptre) == list:
        buildTreeFromXml(sptre[0], tree)
        buildTreeFromXml(sptre[1], tree)
                     
    return tree

def buildRecTreeFromXml(sptre, tree, recType):

    if type(sptre) ==  dict:

        if  len(sptre.keys()) == 1:
            k = list(sptre.keys())[0]
            v = sptre[k]
            if len(v) == 3:
                tree.add_child(name=v['name'])
                clade = tree.search_nodes(name=v["name"])[0]
                event_recs = list(v["eventsRec"].keys())
                event_rec = event_recs[0]
                if recType == "geneSpecie":
                    species_location =v["eventsRec"][event_rec]["speciesLocation"]
                    clade.add_feature("event_rec", event_rec)
                    clade.add_feature("location", species_location)
                elif recType == "transcriptGene":
                    species_location =v["eventsRec"][event_rec]["genesLocation"]
                    clade.add_feature("event_rec", event_rec)
                    clade.add_feature("location", species_location)
                                                    
                return buildRecTreeFromXml(v["clade"], clade, recType)
            else:
                return tree
                
        elif  len(sptre.keys()) == 2:
                tree.add_child(name=sptre['name'])
                clade = tree.search_nodes(name=sptre["name"])[0]
                event_recs = list(sptre["eventsRec"].keys())
                event_rec = event_recs[0]
                if recType == "geneSpecie":                
                    species_location =sptre["eventsRec"][event_rec]["speciesLocation"]
                    clade.add_feature("event_rec", event_rec)
                    clade.add_feature("location", species_location)                
                elif recType == "transcriptGene":                    
                    species_location =sptre["eventsRec"][event_rec]["genesLocation"]
                    clade.add_feature("event_rec", event_rec)
                    clade.add_feature("location", species_location)                                
                return clade
                                
        elif  len(sptre.keys()) == 3:

                tree.add_child(name=sptre['name'])
                clade = tree.search_nodes(name=sptre["name"])[0]
                event_recs = list(sptre["eventsRec"].keys())
                event_rec = event_recs[0]
                if recType == "geneSpecie":                
                    species_location =sptre["eventsRec"][event_rec]["speciesLocation"]
                    clade.add_feature("event_rec", event_rec)
                    clade.add_feature("location", species_location)
                elif recType == "transcriptGene":
                    species_location =sptre["eventsRec"][event_rec]["genesLocation"]
                    clade.add_feature("event_rec", event_rec)
                    clade.add_feature("location", species_location)
                                    
                return buildRecTreeFromXml(sptre["clade"], clade, recType)        
    elif type(sptre) == list:
        buildRecTreeFromXml(sptre[0], tree, recType)
        buildRecTreeFromXml(sptre[1], tree, recType)
                     
    return tree
                            
def parseDictTree(dict_tree):
    if len(dict_tree.keys()) == 1 and list(dict_tree.keys())[0] == "recPhylo":
        recPhylo = dict_tree["recPhylo"]
        if len(recPhylo.keys())==2 and ("spTree" in recPhylo.keys()) and ("recGeneTree" in recPhylo.keys()):
            return recPhylo
        elif len(recPhylo.keys())==2 and ("gnTree" in recPhylo.keys()) and ("recTransTree" in recPhylo.keys()):
            return recPhylo
        else:
            return False
    else:
        return False
        

def makeMapping(rec_gene_tree_nw):
    speciesMapping = {}
    node_mapping_to_parent = {}
    for n in rec_gene_tree_nw.traverse():
        if hasattr(n,"location"):
            node_mapping_to_parent[n.name] = n.get_ancestors()[0].name
            if n.location in speciesMapping.keys():
                speciesMapping[n.location].append([n.name, n.event_rec])            
            else:
                speciesMapping[n.location] = [[n.name, n.event_rec]]
    return speciesMapping, node_mapping_to_parent

def browseSpTree(sp_tree, node_mapping_to_parent):
    for n in sp_tree.traverse():
        node_mapping_to_parent[n.name] = n.get_ancestors()[0].name        
    return node_mapping_to_parent
    
def add_branch_length(tree):
    D = 100.0
    numberOfLeaves = len(tree)
    pas = D/numberOfLeaves
    
    for n in tree.traverse("preorder"):
        if n.is_root():
            n.remainingDistance = D
        elif n.is_leaf():
            ancestor = n.get_ancestors()[0]
            n.dist = ancestor.remainingDistance
        else:
            ancestor = n.get_ancestors()[0]
            ancestorDist = ancestor.remainingDistance
            n.dist = pas
            n.remainingDistance = ancestorDist - pas

                
def parse_rec(my_xml, recType):
    dict_tree = make_dict_from_tree(xml.etree.ElementTree.fromstring(my_xml), recType)
    recPhylo = parseDictTree(dict_tree)

    if recPhylo!= False:
        if recType == "geneSpecie":        
            sptre = recPhylo["spTree"]["phylogeny"]
            recGeneTree = recPhylo["recGeneTree"]["phylogeny"]
        elif recType == "transcriptGene":
            sptre = recPhylo["gnTree"]["phylogeny"]
            recGeneTree = recPhylo["recTransTree"]["phylogeny"]        
        sp_tree = Tree()
        sp_tree_tmp = Tree()       
        rec_gene_tree = Tree()
        sptre_nw = buildTreeFromXml(sptre, sp_tree)
        sp_tree_tmp = sp_tree        
        add_branch_length(sp_tree)
        sp_tree = sp_tree_tmp.write(format=1)
        rec_gene_tree_nw = buildRecTreeFromXml(recGeneTree, rec_gene_tree, recType)    
        speciesMapping, node_mapping_to_parent = makeMapping(rec_gene_tree_nw)
        browseSpTree(sptre_nw, node_mapping_to_parent)
        rec_gene_tree = rec_gene_tree_nw.write(format=1)

        return sp_tree, rec_gene_tree, speciesMapping, node_mapping_to_parent



