#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compute_orthologs.py`` **module description**:

This module computes orthologs groups based on spliced alignments

moduleauthor:: Safa Jammali, 
Université de Sherbrooke Canada

2017-2018
 """

from compute_alignments import *
##############################
### ORTHOLOGY COMPUTATION ####
##############################

def computeOrthology(sourcedata,targetdata,comparisonresults):

    """
     Parameters
     ----------
     sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    comparisonresults: list
	blocs result from used alignment method 




     Returns
     -------
     grouplist:list 
     """
    #print "Computing orthology groups"
    grouplist = []
    
    orthologymatrix = compute_orthology_matrix(sourcedata,targetdata,comparisonresults)
    for gene1 in range(len(orthologymatrix)):
        assigned_to_group = False
        for i in range(len(grouplist)):
            group = grouplist[i]
            group_assigned = False
            for gene2 in group:
                if(orthologymatrix[gene1][gene2] == 1):
                    group_assigned = True
            if(group_assigned):
                assigned_to_group = True
                grouplist[i].append(gene1)
    
        if(not assigned_to_group):
            grouplist.append([gene1])

    grouplist = merge_orthology_groups(grouplist)
    return grouplist

def compute_orthology_matrix(sourcedata,targetdata,comparisonresults):
    """


    Parameters
    ----------
    sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    comparisonresults: list
	blocs result from used alignment method 


    Returns
    -------
    orthologymatrix:matrix


    """
    orthologymatrix = []
    for cds in sourcedata:
        orthologymatrix.append([])
        for cds in sourcedata:
            orthologymatrix[-1].append(0)

    i = 0
    for gene in targetdata:
        geneid,geneseq = gene
        j = 0
        for cds in sourcedata:
            cdsid,cdsseq,cdsgeneid, null = cds
            status, blocklist,splicing_sites,targetcds, texisting = comparisonresults[i][j]
            if(status == STATUS_EXISTING_PROTEIN):
                orthologymatrix[j][texisting] = 1
                orthologymatrix[texisting][j] = 1
            j += 1
        i += 1
    return orthologymatrix

def merge_orthology_groups(grouplist):
    """


     Parameters
     ----------
     grouplist:

     Returns
     -------
     grouplist:
     """
    groups_to_delete = []
    for i in range(len(grouplist)):
        for j in range(i+1,len(grouplist)):
            if(len(list(set(grouplist[i]) & set(grouplist[j]))) != 0):
                grouplist[i] += grouplist[j]
                groups_to_delete.append(grouplist[j])

    newgrouplist=[]
    for sublist in groups_to_delete:
        if sublist not in newgrouplist:
            newgrouplist.append(sublist)
    
    for group in newgrouplist:
        grouplist.remove(group)
    
    return grouplist


def completeOrthology(sourcedata,targetdata,comparisonresults,orthologygroups):
    """


     Parameters
     ----------

    sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    comparisonresults: list
	blocs result from used alignment method 

     orthologygroups:


     Returns
     -------
     sourcedata:
     orthologygroups:

     """
    #print "Adding predicted protein-coding sequence"
    cds2id = []
    cds2geneid = []
    geneidlist = []
    
    nbcdsinitial = len(sourcedata)
    for cds in sourcedata:
        cdsid,cdsseq,cdsgeneid, null = cds
        cds2id.append(cdsid)
        cds2geneid.append(cdsgeneid)

    for k in range(len(orthologygroups)):
        group2geneid = []
        for cds in orthologygroups[k]:
            group2geneid.append(cds2geneid[cds])
        
        for i in range(len(targetdata)):
            gene = targetdata[i]
            geneid,geneseq = gene
            if(geneid not in group2geneid):
                #predicted = False
                for cds in orthologygroups[k]:
                    if(cds < nbcdsinitial):
                        cdsid,null,null,null = sourcedata[cds]
                        status, blocklist,splicing_sites,targetcds, texisting = comparisonresults[i][cds]
                        all_valid_sites = all_valid_splicing_sites(splicing_sites)
                        if((status == STATUS_PREDICTED_PROTEIN or status == STATUS_COMPLETE_CDS) and all_valid_sites):
                            prediction = ""
                            if(status == STATUS_PREDICTED_PROTEIN):
                                prediction = "Protein"
                            else:
                                prediction = "CDS"
                                
                            pid =  "Predicted_"+prediction+"_Group"+str(k)+"_"+cdsid+"_to_"+geneid
                            pseq = targetcds
                            pgeneid = geneid
                            new_prediction = True
                            j = len(sourcedata)-1
                            while(j >= 0 and sourcedata[j][2] == pgeneid):
                                new_prediction = new_prediction and (pseq != sourcedata[j][1])
                                j -= 1
                            if(new_prediction):
                                sourcedata.append([pid,pseq,pgeneid,[]])#TO COMPUTE
                                orthologygroups[k].append(len(sourcedata)-1)
                                
    return sourcedata,orthologygroups


