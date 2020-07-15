#!/usr/bin/python
#-*- coding: utf-8 -*-

import sys
def compute_multi_ortholog(cdsexon,alnfilename, orthofilename):
    alnfile = open(alnfilename, "r")
    mblocklist = []
    idlist = []
    lines = alnfile.readlines()
    i = 0
    while(i < len(lines)):
        line = lines[i]
        if(line[0] == '>'):
            idlist = []
            mblocklist.append({})
            i += 1
            while(len(lines[i]) > 3):
                line = lines[i]
                tab = line.split("\n")[0].split(":")
                seqid = tab[0]
                for k in range(1,len(tab)-1):
                    seqid += ":"+tab[k]
                seqpos = [int(x) for x in tab[-1].split("-")]
                if(seqpos != [0,0]):
                    mblocklist[-1][seqid] = seqpos
                #if("T0" in seqid):
                #    idlist.append(seqid)
                i += 1
        else:
            i += 1

    alnfile.close()

    multigrouplist = []

    idlist = list(cdsexon.keys())
    tmpidlist = []
    for id in idlist:
        tmpidlist.append(id)
            
 
    while(len(tmpidlist) > 0):
        delete = [0]
        id1 = tmpidlist[0]
        group = [id1]
        for j in range(1, len(tmpidlist)):
            id2 = tmpidlist[j]
            ortholog = True
            for mblock in mblocklist:
                if ((id1 in mblock.keys()) != (id2 in mblock.keys())):
                    ortholog = False
                elif(id1 in mblock.keys() and (mblock[id1][1]-mblock[id1][0]-mblock[id2][1]+mblock[id2][0]) % 3 != 0):
                    ortholog = False
            if(ortholog):
                group.append(id2)
                delete.append(j)

        for i in range(len(delete)-1,-1,-1):
            pos = delete[i]
            tmpidlist.remove(tmpidlist[pos])

        multigrouplist.append(group)

                                                        
    multiorthofile = open(orthofilename+"_multi.txt", "w")

    i = 0
    for mgroup in multigrouplist:
        multiorthofile.write(">"+str(i)+"\n")
        for seqid in mgroup:
            multiorthofile.write(seqid+"\n")
        multiorthofile.write("\n")
        i += 1

    multiorthofile.close()

