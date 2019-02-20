# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from Bio import Phylo
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import random
#from cStringIO import StringIO
from io import StringIO
from dash.dependencies import Input, Output
from parseRec import *

app = dash.Dash()

def get_leaf_names(tree):
    leafNames = []
    leaves = tree.get_terminals()
    for leaf in leaves:
        leafNames.append(leaf.name)
    return leafNames
    
def label_tree_internal_nodes(tree):
    i = 999990
    for t in tree.find_clades():
        if t.confidence!= None:
            t.name = t.confidence
            
        if t.name == None:
            t.name = str(i)
            i = i + 1
        t.name = str(t.name)

def get_x_coordinates(tree, width):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
        
    for clade in xcoords.keys():
        if len(clade.clades) ==2:
            xcoords[clade] = [xcoords[clade], xcoords[clade], xcoords[clade], xcoords[clade], xcoords[clade]+width, xcoords[clade]+ width, xcoords[clade]+ width, xcoords[clade]+ width, xcoords[clade]+ width]
        else:
            xcoords[clade] = [xcoords[clade], xcoords[clade]]
			
    return xcoords
    
 
def get_y_coordinates_old(tree, xcoords, width, dist=1.3):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        x1 = xcoords[clade.clades[-1]][0]
        y1 = ycoords[clade.clades[-1]]
        x2 = xcoords[clade.clades[0]][0]
        y2 = ycoords[clade.clades[0]]

        ycoords[clade] = (x1*y2 + x2*y1)/(x1 + x2)
        

    if tree.root.clades:
        calc_row(tree.root)
        
                                
    for clade in tree.find_clades():
        if len(clade.clades) == 2:
            child_bottom = clade.clades[0]
            child_top = clade.clades[-1]
            y_child_bottom = ycoords[child_bottom]
            y_child_top = ycoords[child_top]        
            
            ycoords[clade] = [y_child_bottom, ycoords[clade], ycoords[clade] + width, y_child_top + width, y_child_bottom + width, ycoords[clade], (ycoords[clade] + ycoords[clade] + width)/2.0, ycoords[clade] + width, y_child_top]        
        else: 
            ycoords[clade] = [ycoords[clade], ycoords[clade] + width]        
    
    for clade in tree.find_clades():
        if len(clade.clades) == 2:            
            child_bottom = clade.clades[0]
            child_top = clade.clades[-1]
           
            if len(xcoords[child_bottom])>2:
                xc1 = xcoords[child_bottom][1]
                yc1 = ycoords[child_bottom][1]
                xc2 = xcoords[child_bottom][2]
                yc2 = ycoords[child_bottom][2]
            else:
                xc1 = xcoords[child_bottom][0]
                yc1 = ycoords[child_bottom][0]
                xc2 = xcoords[child_bottom][1]
                yc2 = ycoords[child_bottom][1]
            x1 = xcoords[clade][1]
            y1 = ycoords[clade][1]
            x2 = xcoords[clade][6]
            y2 = ycoords[clade][6]

            yc2 = y2 - ((xc2 -x2)*(y1 - yc1))/(xc1-x1)

            if len(xcoords[child_bottom])>2:
                ycoords[child_bottom][2] =  yc2
            else:
                ycoords[child_bottom][1] =  yc2



            if len(xcoords[child_top])>2:
                xc1 = xcoords[child_top][2]
                yc1 = ycoords[child_top][2]
                xc2 = xcoords[child_top][1]
                yc2 = ycoords[child_top][1]
            else:
                xc1 = xcoords[child_top][1]
                yc1 = ycoords[child_top][1]
                xc2 = xcoords[child_top][0]
                yc2 = ycoords[child_top][0]
            x1 = xcoords[clade][2]
            y1 = ycoords[clade][2]
            x2 = xcoords[clade][6]
            y2 = ycoords[clade][6]
            #print (child_bottom, x2, y2, xc2, yc2, x1, y1, xc1, yc1)
            yc2 = y2 - ((xc2 -x2)*(y1 - yc1))/(xc1-x1)
            #print (yc2)
            #print(ycoords[child_top])
            if len(xcoords[child_top])>2:
                ycoords[child_top][1] =  yc2
            else:
                ycoords[child_top][0] =  yc2
            #print(ycoords[child_top])
            #print ("----------")

    #print (ycoords)     
    return ycoords                 

    

def get_y_coordinates(tree, width, dist):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
  
    for clade in tree.find_clades():
        if len(clade.clades) == 2:
            child_bottom = clade.clades[0]
            child_top = clade.clades[-1]
            y_child_bottom = ycoords[child_bottom]
            y_child_top = ycoords[child_top]        	        	
            ycoords[clade] = [y_child_bottom, ycoords[clade], ycoords[clade] + width, y_child_top + width, y_child_bottom + width, ycoords[clade], (ycoords[clade] + ycoords[clade] + width)/2.0, ycoords[clade] + width, y_child_top]        
        else: 
            ycoords[clade] = [ycoords[clade], ycoords[clade] + width]        
            
    return ycoords
    
    
def get_clade_lines(clade, line_shapes, width, child, line_color='rgb(25,25,25)', line_width=3, x_coords=[0,0], y_coords=[0,0]):
    """define a shape of type 'line', for branch
    """

    x_parent = x_coords[clade]
    x_child = x_coords[child]
    y_parent = y_coords[clade]
    y_child = y_coords[child]	
    
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
                       
    if len(x_coords[clade]) ==  len(x_coords[child]) == 2:
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_parent[0],
                 x1=x_child[0],
                 y1=y_child[0]
                 )
        )
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[1],
                 y0=y_parent[1],
                 x1=x_child[1],
                 y1=y_child[1]
                 )
        )       
    elif len(x_coords[clade]) == 2 and len(x_coords[child]) == 9:                    
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_parent[0],
                 x1=x_child[1],
                 y1=y_child[1]
                 )
        )
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[1],
                 y0=y_parent[1],
                 x1=x_child[2],
                 y1=y_child[2]
                 )
        )  
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[0],
                 y0=y_child[0],
                 x1=x_child[1],
                 y1=y_child[1]
                 )
        ) 
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[2],
                 y0=y_child[2],
                 x1=x_child[3],
                 y1=y_child[3]
                 )
        )             
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[4],
                 y0=y_child[4],
                 x1=x_child[8],
                 y1=y_child[8]
                 )
        )        
    elif len(x_coords[clade]) == 9 and len(x_coords[child]) == 2:                    
        if (y_child[0] > y_parent[0]):        
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[8],
                     y0=y_parent[8],
                     x1=x_child[0],
                     y1=y_child[0]
                     )
            )
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[3],
                     y0=y_parent[3],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
        )
        else:
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[0],
                     y0=y_parent[0],
                     x1=x_child[0],
                     y1=y_child[0]
                     )
            )
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[4],
                     y0=y_parent[4],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            ) 
         
  	    
    elif len(x_coords[clade]) == 9 and len(x_coords[child]) == 9: 
        if (y_child[0] < y_parent[0]):                                 
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[0],
                     y0=y_parent[0],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            )
            
            line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[4],
                 y0=y_parent[4],
                 x1=x_child[2],
                 y1=y_child[2]
                 )
            ) 
        else:
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[3],
                     y0=y_parent[3],
                     x1=x_child[2],
                     y1=y_child[2]
                     )                     
            )
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[8],
                     y0=y_parent[8],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            )            
 
        line_shapes.append(
            dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[8],
                 y0=y_child[8],
                 x1=x_child[4],
                 y1=y_child[4]
                 )
            ) 
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[2],
                 y0=y_child[2],
                 x1=x_child[3],
                 y1=y_child[3]
                 )
        )        
        line_shapes.append(
			dict(type='line',
				layer='below',
				line=dict(color=line_color, width=line_width),
				x0=x_child[0],
				y0=y_child[0],
				x1=x_child[1],
				y1=y_child[1]
				)
			) 
              
        
    return line_shapes


def get_clade_lines_slanted(clade, line_shapes, width, child, line_color='rgb(25,25,25)', line_width=3, x_coords=[0,0], y_coords=[0,0]):
    """define a shape of type 'line', for branch
    """
   
    x_parent = x_coords[clade]
    x_child = x_coords[child]
    y_parent = y_coords[clade]
    y_child = y_coords[child]	
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
              
    if len(x_coords[clade]) ==  len(x_coords[child]) == 2:
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_parent[0],
                 x1=x_child[0],
                 y1=y_child[0]
                 )
        )
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[1],
                 y0=y_parent[1],
                 x1=x_child[1],
                 y1=y_child[1]
                 )
        )       
    elif len(x_coords[clade]) == 2 and len(x_coords[child]) == 9:                    
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_parent[0],
                 x1=x_child[1],
                 y1=y_child[1]
                 )
        )
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[1],
                 y0=y_parent[1],
                 x1=x_child[2],
                 y1=y_child[2]
                 )
        )  
    
    elif len(x_coords[clade]) == 9 and len(x_coords[child]) == 2:                    
        if (y_child[0] > y_parent[0]):        
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[2],
                     y0=y_parent[2],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            )
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[6],
                     y0=y_parent[6],
                     x1=x_child[0],
                     y1=y_child[0]
                     )
        )
        else:
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[1],
                     y0=y_parent[1],
                     x1=x_child[0],
                     y1=y_child[0]
                     )
            )
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[6],
                     y0=y_parent[6],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            ) 
         
  	    
    elif len(x_coords[clade]) == 9 and len(x_coords[child]) == 9: 
        if (y_child[0] < y_parent[0]):                                 
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[1],
                     y0=y_parent[1],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            )
            
            line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[6],
                 y0=y_parent[6],
                 x1=x_child[2],
                 y1=y_child[2]
                 )
            ) 
        else:
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[2],
                     y0=y_parent[2],
                     x1=x_child[2],
                     y1=y_child[2]
                     )                     
            )
            line_shapes.append(
        	    dict(type='line',
                     layer='below',
                     line=dict(color=line_color, width=line_width),
                     x0=x_parent[6],
                     y0=y_parent[6],
                     x1=x_child[1],
                     y1=y_child[1]
                     )
            )            
    else:
        print ("here")
    return line_shapes


def get_clade_lines_slanted_inter(recType,clade, line_shapes, width, child, line_color='rgb(25,25,25)', line_width=3, x_coords=[0,0], y_coords=[0,0], species_mapping={}, node_mapping_to_parent={}):
    """define a shape of type 'line', for branch
    """
   
    x_parent = x_coords[clade]
    x_child = x_coords[child]
    y_parent = y_coords[clade]
    y_child = y_coords[child]	
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
   
    if (len(x_coords[clade]) == len(y_coords[clade]) == 2) and (len(x_coords[child]) == len(y_coords[child]) == 1):
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=(x_parent[0] + x_parent[1])/2.0,
                 y0=(y_parent[0] + y_parent[1])/2.0,
                 x1=x_child[0],
                 y1=y_child[0]
                 )
        )
        
       
    elif len(x_coords[clade]) == len(y_coords[clade]) == len(x_coords[child]) == len(y_coords[child]) == 1:        
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_parent[0],
                 x1=x_parent[0],
                 y1=y_child[0]
                 )
        )
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_child[0],
                 x1=x_child[0],
                 y1=y_child[0]
                 )
        ) 
        
        if species_mapping[child.name][0][1] == "loss" :
            x0=x_child[0],
            x0 = x0[0]
            y0=y_child[0]
            line_shapes.append(
                dict(
                type = 'path',
                path= 'M' + str(x0 - 0.75) + ' ' + str(y0 - 0.75) + 'L' + str(x0+0.75) + ' ' + str(y0+0.75) +'Z',
                fillcolor = 'white',
                line = dict(
                   color = 'red',
                 )            
                )        
            )               
            line_shapes.append(
                dict(
                type = 'path',
                path= 'M' + str(x0+0.75) + ' ' + str(y0-0.75) + 'L' + str(x0 - 0.75) + ' ' + str(y0+0.75) + 'Z',
                fillcolor = 'white',
                line = dict(
                   color = 'red',
                 )            
                )        
            )               
            

    elif (len(x_coords[clade]) == len(y_coords[clade]) == 1) and (len(x_coords[child]) == len(y_coords[child]) == 3) and ((species_mapping[child.name][0][1]=="duplication")or(species_mapping[child.name][0][1]=="creation")):
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_parent[0],
                 x1=x_parent[0],
                 y1=y_child[0]
                 )
        )
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_parent[0],
                 y0=y_child[0],
                 x1=x_child[0],
                 y1=y_child[0]
                 )
        )  
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[0],
                 y0=y_child[0],
                 x1=x_child[0],
                 y1=y_child[1]
                 )
        )          
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[0],
                 y0=y_child[0],
                 x1=x_child[0],
                 y1=y_child[2]
                 )
        )          
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[0],
                 y0=y_child[1],
                 x1=x_child[1],
                 y1=y_child[1]
                 )
        )          
        line_shapes.append(
        	dict(type='line',
                 layer='below',
                 line=dict(color=line_color, width=line_width),
                 x0=x_child[0],
                 y0=y_child[2],
                 x1=x_child[2],
                 y1=y_child[2]
                 )
        )   
        if recType=="geneSpecie":
            x0=x_child[0],
            x0 = x0[0]
            y0=y_child[0]
            line_shapes.append(
                dict(
                type = 'path',
                path= 'M' + str(x0 - 0.5) + ' ' + str(y0 - 0.5) + 'L' + str(x0 - 0.5) + ' ' + str(y0+0.5) + 'L' + str(x0+0.5) + ' ' + str(y0+0.5) + 'L' + str(x0+0.5) + ' ' + str(y0-0.5) + 'Z',
                fillcolor = 'red',
                line = dict(
                   color = 'red',
                 )            
                )        
            )                       
        elif recType=="transcriptGene":             
            x0=x_child[0],
            x0 = x0[0]
            y0=y_child[0]
            line_shapes.append(
                dict(
                type = 'path',
                path= 'M' + str(x0-1.5) + ' ' + str(y0) + 'L' +str(x0) + ' ' + str(y0 - 0.5) + 'L' + str(x0) + ' ' + str(y0+0.5)  +  'Z',
                fillcolor = 'black',
                line = dict(
                   color = 'black',
                 )            
                )        
            )                                  
               
    else:
        print ("here")
    return line_shapes




def draw_clade(clade, slanted, line_shapes, width, child, line_color='rgb(15,15,15)', line_width=3, x_coords=[0,0], y_coords=[0,0]):
    """Recursively draw the tree branches, down from the given clade"""   
    
    if child == None:
        if len(clade.clades) > 0:
            for next_child in clade.clades:
                draw_clade(clade, slanted, line_shapes, width, next_child, line_color, line_width, x_coords, y_coords)    
    else:
        if slanted == False:
            line_shapes = get_clade_lines(clade, line_shapes, width, child, line_color, line_width, x_coords, y_coords)
        else:
            line_shapes = get_clade_lines_slanted(clade, line_shapes, width, child, line_color, line_width, x_coords, y_coords)
            
        if len(child.clades) > 0:
            for next_child in child.clades:
                draw_clade(child, slanted, line_shapes, width, next_child, line_color, line_width, x_coords, y_coords)    


def draw_clade_inter(recType, clade, slanted, line_shapes, width, child, line_color='rgb(255,0,0)', line_width=3, x_coords=[0,0], y_coords=[0,0], species_mapping={}, node_mapping_to_parent={}):
    """Recursively draw the tree branches, down from the given clade"""   
    
    if child == None:
        if len(clade.clades) > 0:
            for next_child in clade.clades:
                draw_clade_inter(recType,clade, slanted, line_shapes, width, next_child, line_color, line_width, x_coords, y_coords, species_mapping, node_mapping_to_parent)    
    else:
        if slanted == False:
            line_shapes = get_clade_lines_slanted_inter(recType,clade, line_shapes, width, child, line_color, line_width, x_coords, y_coords, species_mapping, node_mapping_to_parent)
        else:
            line_shapes = get_clade_lines_slanted_inter(recType, clade, line_shapes, width, child, line_color, line_width, x_coords, y_coords, species_mapping, node_mapping_to_parent)
            
        if len(child.clades) > 0:
            for next_child in child.clades:
                draw_clade_inter(recType, child, slanted, line_shapes, width, next_child, line_color, line_width, x_coords, y_coords, species_mapping, node_mapping_to_parent)    

            
def random_color():
    # Make a list of colors to picvk from
    colors = ["red", "green", "blue", "orange", "purple", "pink", "yellow", "black", "gray", "sliver", "violet",
              "yellowgreen", "turquoise", "sienna", "salmon"]
    color = random.choice(colors)
    return tuple(color)
    
def addInternalTreeTransGene(clade, rec_gene_tree, species_mapping, x_coords, y_coords, x_coords_internal, y_coords_internal, x_name_coords, width, parent, tree, node_mapping_to_parent, nb_inter_elt ={}):        
    if clade.name in species_mapping.keys():

        for inter_clade in species_mapping[clade.name]:            
            nb_inter_nodes = len(species_mapping[clade.name])     
            
            if len(y_coords[clade]) == 2:
                a_y = y_coords[clade][0]
                b_y = y_coords[clade][1]
                
                a_x = x_coords[clade][0]
                b_x = x_coords[clade][1]
                
            elif len(y_coords[clade]) == 9:
                a_y = y_coords[clade][1]
                b_y = y_coords[clade][2]

                a_x = x_coords[clade][0] + width/2.0
                b_x = x_coords[clade][1] + width/2.0
                
            else:
                print ("ERROR")
                exit()
            
            y = []
            x = []

            nb_inter_nodes_dup = 0
            for event in species_mapping[clade.name]:
                if event[1]  == "creation":
                    nb_inter_nodes_dup +=1
            
            nb_inter_elt[clade.name] = [nb_inter_nodes, nb_inter_nodes_dup]

            j= 0
            for i in range(nb_inter_nodes):
                if species_mapping[clade.name][i][1] == "creation":
                    y.append(a_y + (j+1)*(b_y-a_y)/(nb_inter_nodes_dup +1))
                    j +=1
                else:
                    y.append(a_y + (i+1-j)*(b_y-a_y)/(nb_inter_nodes - nb_inter_nodes_dup +1))
                 
                x.append(a_x + (i+1)*(b_x-a_x)/(nb_inter_nodes+1))
                 
            y_coords_internal[clade] = y

            x_shift = width/(nb_inter_nodes + 1)
                        
            
            if inter_clade[1] == "creation":
                if parent == None:
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(clade)/2.0 + x_shift * (len(x_coords_internal[clade]) + 1))
                    else:
                        x_coords_internal[clade]=  [tree.distance(clade)/2.0 + x_shift]                                    
                else:
                    
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0 + x_shift * (len(x_coords_internal[clade]) + 1))
                    else:
                        x_coords_internal[clade]=  [tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0 + x_shift]                
            
            elif inter_clade[1] == "leaf":
                if clade in x_coords_internal.keys():
                    x_coords_internal[clade].append(x_coords[clade][1])
                else:
                    x_coords_internal[clade] = [x_coords[clade][1]]                                           
            
            elif inter_clade[1] == "duplication":
                if clade in x_coords_internal.keys():
                    x_coords_internal[clade].append(x_coords[clade][1] + x_shift * (len(x_coords_internal[clade]) + 1))
                else:
                    x_coords_internal[clade] = [x_coords[clade][1] + x_shift]                
                
            elif inter_clade[1] == "loss":
                if parent == None:
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(clade)/2.0 + x_shift * (len(x_coords_internal[clade]) + 1))
                    else:
                        x_coords_internal[clade]=  [tree.distance(clade)/2.0 + x_shift]                                                    
                else:
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0)
                    else:
                        x_coords_internal[clade]=  [tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0]
                
    if len(clade.clades)>0:
        for next_child in clade.clades:  
            if clade.branch_length == None:
                clade.branch_length = 0

            addInternalTreeTransGene(next_child, rec_gene_tree, species_mapping, x_coords, y_coords, x_coords_internal, y_coords_internal, x_name_coords, width, clade, tree, node_mapping_to_parent, nb_inter_elt)


def addInternalTreeGeneSpecie(clade, rec_gene_tree, species_mapping, x_coords, y_coords, x_coords_internal, y_coords_internal, x_name_coords, width, parent, tree, node_mapping_to_parent, nb_inter_elt ={}):        
    if clade.name in species_mapping.keys():

        for inter_clade in species_mapping[clade.name]:            
            nb_inter_nodes = len(species_mapping[clade.name])     
            
            if len(y_coords[clade]) == 2:
                a_y = y_coords[clade][0]
                b_y = y_coords[clade][1]
                
                a_x = x_coords[clade][0]
                b_x = x_coords[clade][1]
                
            elif len(y_coords[clade]) == 9:
                a_y = y_coords[clade][1]
                b_y = y_coords[clade][2]

                a_x = x_coords[clade][0] + width/2.0
                b_x = x_coords[clade][1] + width/2.0
                
            else:
                print ("ERROR")
                exit()
            
            y = []
            x = []

            nb_inter_nodes_dup = 0
            for event in species_mapping[clade.name]:
                if event[1]  == "duplication":
                    nb_inter_nodes_dup +=1
            
            nb_inter_elt[clade.name] = [nb_inter_nodes, nb_inter_nodes_dup]

            j= 0
            for i in range(nb_inter_nodes):
                if species_mapping[clade.name][i][1] == "duplication":
                    y.append(a_y + (j+1)*(b_y-a_y)/(nb_inter_nodes_dup +1))
                    j +=1
                else:
                    y.append(a_y + (i+1-j)*(b_y-a_y)/(nb_inter_nodes - nb_inter_nodes_dup +1))
                 
                x.append(a_x + (i+1)*(b_x-a_x)/(nb_inter_nodes+1))
                 
            y_coords_internal[clade] = y

            x_shift = width/(nb_inter_nodes + 1)
                        
            
            if inter_clade[1] == "duplication":
                if parent == None:
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(clade)/2.0 + x_shift * (len(x_coords_internal[clade]) + 1))
                    else:
                        x_coords_internal[clade]=  [tree.distance(clade)/2.0 + x_shift]                                    
                else:
                    
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0 + x_shift * (len(x_coords_internal[clade]) + 1))
                    else:
                        x_coords_internal[clade]=  [tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0 + x_shift]                
            
            elif inter_clade[1] == "leaf":
                if clade in x_coords_internal.keys():
                    x_coords_internal[clade].append(x_coords[clade][1])
                else:
                    x_coords_internal[clade] = [x_coords[clade][1]]                                           
            
            elif inter_clade[1] == "speciation":
                if clade in x_coords_internal.keys():
                    x_coords_internal[clade].append(x_coords[clade][1] + x_shift * (len(x_coords_internal[clade]) + 1))
                else:
                    x_coords_internal[clade] = [x_coords[clade][1] + x_shift]                
                
            elif inter_clade[1] == "loss":
                if parent == None:
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(clade)/2.0 + x_shift * (len(x_coords_internal[clade]) + 1))
                    else:
                        x_coords_internal[clade]=  [tree.distance(clade)/2.0 + x_shift]                                                    
                else:
                    if clade in x_coords_internal.keys():
                        x_coords_internal[clade].append(tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0)
                    else:
                        x_coords_internal[clade]=  [tree.distance(parent) + (tree.distance(clade) - tree.distance(parent))/2.0]
                
    if len(clade.clades)>0:
        for next_child in clade.clades:  
            if clade.branch_length == None:
                clade.branch_length = 0

            addInternalTreeGeneSpecie(next_child, rec_gene_tree, species_mapping, x_coords, y_coords, x_coords_internal, y_coords_internal, x_name_coords, width, clade, tree, node_mapping_to_parent, nb_inter_elt)
                    
                    
def create_tree(recData, slanted, recType):
    dist = 15

    sptree, rec_gene_tree, species_mapping, node_mapping_to_parent = parse_rec(recData, recType)    
    
    tree = read_tree_nw(sptree)    
    
    label_tree_internal_nodes(tree) 

    dic = tree.depths()
    val1= max(dic.items(), key = lambda x: x[1])[1]
    maxheight = tree.count_terminals()  
    dic = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))
    val2= max(dic.items(), key = lambda x: x[1])[1]
    width = (val2*2.0)/val1 + 6.5
    
    x_coords = get_x_coordinates(tree, width)
    #y_coords = get_y_coordinates_old(tree, x_coords, width, dist)
    y_coords = get_y_coordinates(tree, width, dist)              
    
    line_shapes = []
    child = None
    
    draw_clade(tree.root, slanted, line_shapes, 20*width, child, line_color=random_color(), line_width=3, x_coords=x_coords, y_coords=y_coords)
    
    x_coords_internal = {}
    y_coords_internal = {}
    x_name_coords = {}
    
    for clade in x_coords.keys():
        x_name_coords[clade.name] = clade
    
    if recType == "geneSpecie":
        addInternalTreeGeneSpecie(tree.root, rec_gene_tree, species_mapping, x_coords, y_coords, x_coords_internal, y_coords_internal, x_name_coords, width, None, tree, node_mapping_to_parent, nb_inter_elt = {})
    elif recType == "transcriptGene":
        addInternalTreeTransGene(tree.root, rec_gene_tree, species_mapping, x_coords, y_coords, x_coords_internal, y_coords_internal, x_name_coords, width, None, tree, node_mapping_to_parent, nb_inter_elt = {})
        
    child = None
    x_coords_internal[tree.root] = x_coords[tree.root]
    y_coords_internal[tree.root] = y_coords[tree.root]   

    draw_clade_inter(recType,tree.root, slanted, line_shapes, 20*width, child, line_color="red", line_width=3, x_coords=x_coords_internal, y_coords=y_coords_internal, species_mapping =species_mapping, node_mapping_to_parent=node_mapping_to_parent)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []
    nodes = []
    
    intermediate_node_color = 'rgb(100,100,100)'

    label_legend = []

    for k in x_coords:
        if k in x_coords_internal.keys():
            for j in range(len(x_coords_internal[k])):
                X.append(x_coords_internal[k][j])
                Y.append(y_coords_internal[k][j]) 
                label_legend.append(k.name)
             
    for cl in my_tree_clades:
        for i in range(len(x_coords[cl])):
            if len(x_coords[cl]) == 2:
                X.append(x_coords[cl][0])
                Y.append(y_coords[cl][0])            
                text.append(cl.name)
                label_legend.append(cl.name)
                label_legend.append(cl.name)
                X.append(x_coords[cl][1])
                Y.append(y_coords[cl][1])            
                text.append(cl.name)
            else:                	
                if slanted == False:
                    if i== 5 or i==6 or i==7: #False : 
                        pass
                    else:
                        X.append(x_coords[cl][i])
                        Y.append(y_coords[cl][i])            
                        text.append(cl.name)
                        label_legend.append(cl.name)
                else:
                    if i== 5 or i==6 or i==7: #False : 
                        pass
                    else:
                        X.append(x_coords[cl][i])
                        Y.append(y_coords[cl][i])            
                        text.append(cl.name) 
                        label_legend.append(cl.name)   

 

    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
                )

    #label_legend = get_leaf_names(tree)
    color = [intermediate_node_color] * len(X)
    
    """
    print (label_legend, len(label_legend), len(X), len(Y))
    
    for elt in label_legend:
        node = dict(type='scatter',
                x=X,
                y=Y,
                mode='markers',
                marker=dict(color=color, size=5),
                text=elt,  # vignet information of each node
                hoverinfo='hoverinfo',
                name=elt
                )
        nodes.append(node)
    """
        
    nodes = [dict(type='scatter',
            x=X,
            y=Y,
            mode='markers',
            #marker=dict(color=color, size=5),
            text=label_legend,  # vignet information of each node
            hoverinfo='text',
            name=label_legend,
            marker = {
                "size" : 5,
                "color" : 'rgb(100,100,100)',
                "line" : {
                    "width" : 2
                }
            })
            ]

    if recType == "geneSpecie":
        graphTitle = "Gene-Specie Reconciliation"
    elif recType == "transcriptGene":
        graphTitle = "Protein-Gene Reconciliation"
        
    layout = dict(title=graphTitle,
                  #paper_bgcolor='rgba(0,0,0,0)',
                  dragmode="select",
                  font=dict(family='Balto', size=14),
                  #width="100%",
                  height=800,
                  autosize=True,
                  showlegend=True,                               
                  xaxis=dict(
                        range = [0, 100],
						titlefont=dict(
							family='Arial, sans-serif',
							size=18,
							color='lightgrey'
						),
						showticklabels=True,
						tickangle=45,
						tickfont=dict(
							family='Old Standard TT, serif',
							size=14,
							color='black'
						),
						exponentformat='e',
						showexponent='all',	                  
						showline=True,
						zeroline=False,
						showgrid=True,  # To visualize the vertical lines
						ticklen=4,
						#showticklabels=True,
						title='branch length'),
                    yaxis=dict(
                        range = [-50, 20],                    
                        scaleanchor = "x",
                        automargin=True,
                        domain = [1,1],                        
                        title = "",
						showline=False,
						zeroline=False,
						showticklabels=False,						
						showgrid=False,  # To visualize the vertical lines
                        titlefont=dict(
	                        family='Arial, sans-serif',
	                        size=18,
	                        color='lightgrey'
						)),

                    hovermode='closest',
                    shapes=line_shapes,
                    plot_bgcolor='rgb(250,250,250)',
                    paper_bgcolor='rgb(250,250,250)',                  
                    legend={'x': 0, 'y': 1}
                  )

    fig = dict(data=nodes, layout=layout)  
    return fig


def read_treefile(filename):
    tree = Phylo.read(filename, "newick")
    return tree


def read_tree_nw(tree_nw):
    tree = Phylo.read(StringIO(tree_nw), "newick")

    return tree
    
    
def read_metadata(filename):
    df = pd.read_csv(filename)
    return df


def create_title(virus, nb_genome):
    graph_title = "Phylogeny of " + virus.title() + " Virus<br>" + str(
        nb_genome) + " genomes colored according to region and country"
    return graph_title

def random_color():
    # Make a list of colors to picvk from
    colors = ["red", "green", "blue", "orange", "purple", "pink", "yellow", "black", "gray", "sliver", "violet",
              "yellowgreen", "turquoise", "sienna", "salmon"]
    color = random.choice(colors)
    return tuple(color)

def dataFromFile(path):
    file = open(path, "r")
    return file.read()
    

slanted = False

recProteinTree = dataFromFile("./datas/recProteinGene.xml")

recGeneSpecie = dataFromFile("./datas/recGeneSpece.xml")

figProteinGene = create_tree(recProteinTree, slanted, "transcriptGene")

figGeneSpecie = create_tree(recGeneSpecie, slanted, "geneSpecie")



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']               
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
    html.H1(children='Double reconciliation visualiszation tool',style={'textAlign': 'center', 'color': '#000000'}),

    html.H6(children='''
       This web application allows the double representation of a gene tree embeded inside one specie tree and transcript tree embeded inside one gene tree.
    ''', style={'textAlign': 'center', 'color': '#000000'}),

    dcc.Dropdown(
        id='my-dropdown',
        options=[
            {'label': 'Slanted Tree', 'value': 'True'},
            {'label': 'No Slanted', 'value': 'False'},
        ],
        style = {"display": "none"},
        value='slanted'
    ),
    dcc.Textarea(
        id = "proteinGene", 
        title = "Protein Gene Reconciliation",   
        placeholder='Enter a value...',
        value=recProteinTree,
        style={'width': '45%','height':'150', 'margin-right': '9%'}
    ),    
    dcc.Textarea(
        id = 'geneSpecie',    
        title = "Gene Specie Reconciliation",        
        placeholder='Enter a value...',
        value=recGeneSpecie,
        style={'width': '45%','height':'150'}
    ),     


    html.Button('Compute and draw', id='button', style={'textAlign': 'center', 'color': 'green', 'margin-left': '43%'}),

    
    dcc.Graph(
    
        #style={
         #   'width': 600,
        #    'height':600
        #},
        id='figProteinGene',
        #Protein Gene Reconciliation',
        figure=figProteinGene
    ),
    dcc.Graph(
    
        #style={
         #   'width': 600,
        #    'height':600
        #},
        id='figGeneSpecie',
        #Gene Specie Reconciliation
        figure=figGeneSpecie
    )    
])


               
@app.callback(Output('figGeneSpecie', 'figure'),
    [Input('button', 'n_clicks'),
    Input('proteinGene', 'value'),
    Input('geneSpecie', 'value')])
def update_output(n_clicks, input2, input3):
    figGeneSpecie = create_tree(input3, False, "geneSpecie")
    return figGeneSpecie
   

@app.callback(Output('figProteinGene', 'figure'),
    [Input('button', 'n_clicks'),
    Input('proteinGene', 'value'),
    Input('geneSpecie', 'value')])
def update_output(n_clicks, input2, input3):
    figProteinGene = create_tree(input2, False, "transcriptGene")
    return figProteinGene   
    
if __name__ == '__main__':
    app.run_server(debug=True)
    
    