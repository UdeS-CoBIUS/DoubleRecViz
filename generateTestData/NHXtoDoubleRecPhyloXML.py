# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Aida Ouangraoua
@date : July 2020
@location : University of Sherbrooke

'''

import os
from functions import *

def NHXtoDoubleRecPhyloXML(file_nhx, file_species_tree):
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
                        command = "python2.7 recPhyloXML-master-20-07-2020/python/NHXtoRecPhyloXML.py -g "+tmp_file_nhx+ " -s "+ file_species_tree+" -o "+ tmp_file_xml + " --include.species 2>/dev/null >/dev/null"
                        first_genetree = False
                    else:
                        command = "python2.7 recPhyloXML-master-20-07-2020/python/NHXtoRecPhyloXML.py -g "+tmp_file_nhx+ " -s "+ file_species_tree+" -o "+ tmp_file_xml + " 2>/dev/null >/dev/null"
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

                    command = "python2.7 recPhyloXML-master-20-07-2020/python/NHXtoRecPhyloXML.py -g "+tmp_file_nhx+ " -s "+ tmp_genetree_file_nw +" -o "+ tmp_file_xml + " 2>/dev/null >/dev/null"
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
