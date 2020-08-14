# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

'''

import argparse
import sys
import os

Dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(Dir+'/../Utils/')

from utils import *

def NHXtoRecPhyloXML(file_nhx, file_tree,type):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    if(dir_name ==  ""):
        dir_name =  "."
    NHXtoRecPhyloXML = dir_name+"/recPhyloXML-master-20-07-2020/python/NHXtoRecPhyloXML.py"
    first_nhx = True
    tmp_file_nhx = "tmp_file.nhx"
    tmp_file_xml = "tmp_file.xml"
    s_xml = ""
    try :
        f = open(file_nhx, 'r')
        for line in f:
            write_tree(line,tmp_file_nhx, 'w')
            # save gene tree at nhx format with root
            s = remove_content_tree_nhx(tmp_file_nhx)
            if(type == "gene"):
                s = add_root_nhx(s,"G")
            else:
                s = add_root_nhx(s,"T")
            write_tree(s,tmp_file_nhx, 'w')
                    
            if(first_nhx == True):
                command = "python2.7 "+NHXtoRecPhyloXML+" -g "+tmp_file_nhx+ " -s "+ file_tree+" -o "+ tmp_file_xml + " --include.species 2>/dev/null >/dev/null"
                first_nhx = False
            else:
                command = "python2.7 "+NHXtoRecPhyloXML+" -g "+tmp_file_nhx+ " -s "+ file_tree+" -o "+ tmp_file_xml + " 2>/dev/null >/dev/null"
                
            os.system(command)
            
            if(type == "gene"):
                # format recGeneTreeXML
                s = read_recGeneTreeXML(tmp_file_xml)
                s_xml += s
            else:
                # format recTransTreeXML
                s = read_recTransTreeXML(tmp_file_xml)
                s_xml += s
        s_xml = "<recPhylo>\n" + s_xml + "</recPhylo>\n"
        os.system("rm " +tmp_file_nhx)
        os.system("rm " +tmp_file_xml)
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
    return s_xml



def NHXtoDoubleRecPhyloXML(file_nhx, file_species_tree):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    if(dir_name ==  ""):
        dir_name =  "."
    NHXtoRecPhyloXML = dir_name+"/recPhyloXML-master-20-07-2020/python/NHXtoRecPhyloXML.py"
    new_genetree = True
    first_genetree = True
    tmp_file_nhx = "tmp_file.nhx"
    tmp_genetree_file_nhx = "tmp_genetree_file.nhx"
    tmp_genetree_file_nw = "tmp_genetree_file.nw"
    s_xml = ""
    try :
        f = open(file_nhx, 'r')
        for line in f:
            if(line == "//\n"):
                new_genetree = True
            else:
                write_tree(line,tmp_file_nhx, 'w')
                tmp_file_xml = "tmp_file.xml"
                if(new_genetree == True):
                    # save gene tree at nhx format with root
                    s = remove_content_tree_nhx(tmp_file_nhx)
                    s = add_root_nhx(s,"G")
                    write_tree(s,tmp_file_nhx, 'w')
                    os.system("cp "+tmp_file_nhx+" "+tmp_genetree_file_nhx)
                    
                    # save gene tree at nw format with root
                    s = remove_content_tree_nw(tmp_genetree_file_nhx)
                    s = add_root_nw(s,"G")
                    write_tree(s,tmp_genetree_file_nw, 'w')
                    
                    if(first_genetree == True):
                        command = "python2.7 "+NHXtoRecPhyloXML+" -g "+tmp_file_nhx+ " -s "+ file_species_tree+" -o "+ tmp_file_xml + " --include.species 2>/dev/null >/dev/null"
                        first_genetree = False
                    else:
                        command = "python2.7 "+NHXtoRecPhyloXML+" -g "+tmp_file_nhx+ " -s "+ file_species_tree+" -o "+ tmp_file_xml + " 2>/dev/null >/dev/null"
                    os.system(command)
                    
                    # format recGeneTreeXML
                    s = read_recGeneTreeXML(tmp_file_xml)
                    s_xml += s
                    new_genetree = False
                else:
                    # save transcript tree at nhx format with root
                    s = remove_content_tree_nhx(tmp_file_nhx)
                    s = add_root_nhx(s,"T")
                    write_tree(s,tmp_file_nhx, 'w')

                    command = "python2.7 "+NHXtoRecPhyloXML+" -g "+tmp_file_nhx+ " -s "+ tmp_genetree_file_nw +" -o "+ tmp_file_xml + " 2>/dev/null >/dev/null"
                    os.system(command)
                    
                    # format recTransTreeXML
                    s = read_recTransTreeXML(tmp_file_xml)
                    s_xml += s
        s_xml = "<doubleRecPhylo>\n" + s_xml + "</doubleRecPhylo>\n"
        os.system("rm " +tmp_file_nhx)
        os.system("rm " +tmp_genetree_file_nhx)
        os.system("rm " +tmp_genetree_file_nw)
        os.system("rm " +tmp_file_xml)
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
    return s_xml

#####################
### Main ############

def build_arg_parser():
    parser = argparse.ArgumentParser(description="Convert NHX to DoubleRecPhyloXML")
    parser.add_argument('-t', '--type', help="type of reconciled trees: genespecies, transcriptgene, or double (required)" )
    parser.add_argument('-s', '--speciesTree', help="name of the species tree file at Newick format (required if -t double, or -t genespecies)")
    parser.add_argument('-rg', '--recGeneTree', help="name of the file containing reconciled gene trees at NHX format (required if -t genespecies)")
    parser.add_argument('-rgt', '--recGeneTransTree', help="name of the file containing reconciled gene and transcript trees at NHX format (required if -t double)")
    parser.add_argument('-g', '--geneTree', help="name of the gene tree file at Newick format (required if -t transcriptgene)")
    parser.add_argument('-rt', '--recTransTree', help="name of the file containing reconciled transcript trees at NHX format (required if -t trancriptgene)")
    parser.add_argument('-o', '--output', help="name of the output file (optional, default is reconciliation_file + \".xml\")")
    return parser


if __name__ == '__main__':
    parser = build_arg_parser()
    args = parser.parse_args()

    o = args.output
    type = args.type
    if(type == None):
        print("Argument -t <type> is required")
        parser.print_help()
    if(type == "double"):
        s = args.speciesTree
        if(s == None):
            print("Argument -s <speciesTree> is required")
            parser.print_help()
        else:
            rgt = args.recGeneTransTree
            if(rgt == None):
                print("Argument -rgt <recGeneTransTree> is required")
                parser.print_help()
            else:
                if(o == None):
                    o = rgt+".xml"
                b = NHXtoDoubleRecPhyloXML(rgt, s)
                write_tree(b,o,'w')
    elif(type == "genespecies"):
        s = args.speciesTree
        if(s == None):
            print("Argument -s <speciesTree> is required")
            parser.print_help()
        else:
            rg = args.recGeneTree
            if(rg == None):
                print("Argument -rg <recGeneTree> is required")
                parser.print_help()
            else:
                if(o == None):
                    o = rg+".xml"
                b = NHXtoRecPhyloXML(rg, s, "gene")
                write_tree(b,o,'w')
    elif(type == "transcriptgene"):
        g = args.geneTree
        if(g == None):
            print("Argument -g <geneTree> is required")
            parser.print_help()
        else:
            rt = args.recTransTree
            if(rt == None):
                print("Argument -rt <recTransTree> is required")
                parser.print_help()
            else:
                if(o == None):
                    o = rt+".xml"
                b = NHXtoDoubleRecPhyloXML(rt, g,"transcript")
                write_tree(b,o,'w')
    else:
        print("Argument -t <type> : type must be double, gene, or transcript")
        parser.print_help()
