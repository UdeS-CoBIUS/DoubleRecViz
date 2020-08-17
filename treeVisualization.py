# -*- coding: utf-8 -*-
#!/usr/bin/python3.7
'''
@author : Esaie Kuitche
@date : April 2020
@update : Yanchun Qi
@date : July 2020

@location : University of Sherbrooke

'''
import dash
import dash_core_components as dcc
import dash_html_components as html
from Bio import Phylo
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import random
#from cStringIO import StringIO #python2.7
from io import StringIO #python3
from dash.dependencies import Input, Output, State
from parseRec import *
from doubleRecPhylo2recPhylo import *
import base64
import lxml


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
		pass
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
			

	elif (len(x_coords[clade]) == len(y_coords[clade]) == 1) and (len(x_coords[child]) == len(y_coords[child]) == 3) and ((species_mapping[child.name][0][1]=="speciation") or (species_mapping[child.name][0][1]=="duplication")or(species_mapping[child.name][0][1]=="creation")):
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
		pass
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
			
			elif (inter_clade[1] == "duplication") or (inter_clade[1] == "speciation"):
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
					
					
def create_tree(recData, slanted, recType, textcolor):
	# change: fit the window
	dist = 9

	sptree, rec_gene_tree, species_mapping, node_mapping_to_parent = parse_rec(recData, recType)	

	nb_dup = 0
	nb_loss = 0
	nb_creat = 0
	for node, events in species_mapping.items():
		for event in events:
			# print ("Debug 1st", event)
			if event[1] == "duplication":
				nb_dup +=1
			elif event[1] == "loss":
				nb_loss +=1
			elif event[1] == "creation":
				nb_creat +=1

	if recType == "geneSpecie":
		text_legend = "The Reconciliation cost is : " + str(nb_dup) + " duplication(s) + " +  str(nb_loss) + " loss = " + str(nb_loss+nb_dup) + " event(s)"
	elif recType == "transcriptGene":
		text_legend = "The Reconciliation cost is : " + str(nb_creat) + " creation(s) + " +  str(nb_loss) + " loss = " + str(nb_loss+nb_creat) + " event(s)"

	reconciliation_cost = 0
	
	tree = read_tree_nw(sptree)	
	
	label_tree_internal_nodes(tree) 

	dic = tree.depths()
	val1= max(dic.items(), key = lambda x: x[1])[1]
	maxheight = tree.count_terminals()  
	dic = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))
	val2= max(dic.items(), key = lambda x: x[1])[1]
	width = (val2*2.0)/val1 + 6.5
	
	x_coords = get_x_coordinates(tree, width)
	y_coords = get_y_coordinates(tree, width, dist)			  
	
	line_shapes = []
	child = None
	if recType == "geneSpecie":	
		draw_clade(tree.root, slanted, line_shapes, 20*width, child, line_color="black", line_width=3, x_coords=x_coords, y_coords=y_coords)
	elif recType == "transcriptGene":
		draw_clade(tree.root, slanted, line_shapes, 20*width, child, line_color="red", line_width=3, x_coords=x_coords, y_coords=y_coords)	        

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
	
	draw_clade_inter(recType,tree.root, slanted, line_shapes, 20*width, child, line_color=textcolor, line_width=3, x_coords=x_coords_internal, y_coords=y_coords_internal, species_mapping =species_mapping, node_mapping_to_parent=node_mapping_to_parent)
	my_tree_clades = x_coords.keys()
	X = []
	Y = []
	text = []
	nodes = []
	
	intermediate_node_color = 'rgb(100,100,100)'

	label_legend = []
	

	rec_gene_tree

	rec_tree = read_tree_nw(rec_gene_tree)	

	terminal_gene = rec_tree.get_terminals()

	terminalNodeNamePara = []
	for e in terminal_gene:
		terminalNodeNamePara.append(e.name)

	
	terminal = tree.get_terminals()
	terminalNodeNameHost = []
	for e in terminal:
		terminalNodeNameHost.append(e.name)
	# print(terminalNodeNameHost, "______________________________")
	# print(x_coords)
	# print(x_coords_internal.keys())
	"""	
	for node in tree.get_descendants():
		# Do some analysis on node
		terminalNodeNameHost.append(node.name)
	"""
	for k in x_coords:		
		if k in x_coords_internal.keys() and k.name != "999990":

			events = species_mapping[k.name]
			selected_events_name = []
			index_selected = []
			for e in events:
				# print("Debug 2nd: ", k.name, "\t",e[0], e[1])
				if e[1] == "leaf" or e[1] == "loss" or e[1] == "duplication" or e[1] == "creation" or e[1] == "speciation": 
					selected_events_name.append(e[0])
					index_selected.append(events.index(e))
				#else:
				#	selected_events_name.append(e[0])
				#	index_selected.append(events.index(e))  

			if len(x_coords_internal[k])==1:
					X.append(x_coords_internal[k][0])
					Y.append(y_coords_internal[k][0]) 
					leaf_name = selected_events_name.pop()
					label_legend.append("\t" + "\t" + leaf_name)
					
			elif len(x_coords_internal[k])==3:
					X.append(x_coords_internal[k][2])
					Y.append(y_coords_internal[k][2]) 
					leaf_name = selected_events_name.pop()
					label_legend.append("\t" + "\t" + leaf_name)

					X.append(x_coords_internal[k][1])
					Y.append(y_coords_internal[k][1]) 
					leaf_name = selected_events_name.pop()
					label_legend.append("\t" + "\t" + leaf_name)

					X.append(x_coords_internal[k][0])
					Y.append(y_coords_internal[k][0]) 
					leaf_name = selected_events_name.pop()
					label_legend.append("\t" + "\t" + leaf_name)

	"""    			    
			for j in range(len(x_coords_internal[k])):
				if j in index_selected:
					print("________________", j, "______________", events[j])
					if k.name not in terminalNodeNameHost:
						X.append(x_coords_internal[k][j])
						Y.append(y_coords_internal[k][j]) 
						leaf_name = selected_events_name.pop()
						label_legend.append("\t" + "\t" + leaf_name)
						print("++++++++++++",leaf_name)					
						#pass
					else:				
						X.append(x_coords_internal[k][j])
						Y.append(y_coords_internal[k][j]) 
						leaf_name = selected_events_name.pop()
						label_legend.append("\t" + "\t" + leaf_name)
						print("_______",leaf_name)
	"""
	host_leaves_names = []
	for x in terminalNodeNameHost:
		host_leaves_names.append({"label": x, 'value': x})					
	
	parasite_leaves_names = []
	for x in terminalNodeNamePara:
		
		if x in ["LOSS", "Loss", "loss"]:
			continue
		else:
			parasite_leaves_names.append({"label": x, 'value': x})	


	for cl in my_tree_clades:
		if cl.name == "999990":
			pass
		else:
			if len(x_coords[cl]) == 2:
				label_legend.append(cl.name)
				X.append(x_coords[cl][1])
				Y.append(y_coords[cl][1])			
				text.append(cl.name)
			elif len(x_coords[cl]) > 2:
				label_legend.append(cl.name)
				X.append(x_coords[cl][1])
				Y.append(y_coords[cl][1])			
				text.append(cl.name)		        
				"""
				label_legend.append(cl.name)
				X.append(x_coords[cl][2])
				Y.append(y_coords[cl][2])			
				text.append(cl.name)
				
				label_legend.append(cl.name)
				X.append(x_coords[cl][3])
				Y.append(y_coords[cl][3])			
				text.append(cl.name)						        
				"""
				
	axis = dict(showline=False,
				zeroline=False,
				showgrid=False,
				showticklabels=False,
				title=''  # y title
				)


	color = [intermediate_node_color] * len(X)
	
	
	nodes = [dict(type='scatter',
			x=X,
			y=Y,
			mode='markers+text',
			#marker=dict(color=color, size=5),
			text=label_legend,  # vignet information of each node
			hoverinfo='text',			
			textposition='middle right',
			#font=dict(family='Balto', size=20, color=color),
			# name=label_legend,
			textfont = {
				"size" : 14,
				"color" : "black"
			},			
			marker = {
				"size" : 1,
				"color" : 'rgb(100,100,100)',
				"line" : {
					"width" : 1
				}
			}
			)
			]

	#nodes = []

	if recType == "geneSpecie":
		graphTitle = "Gene-species reconciliation"
	elif recType == "transcriptGene":
		graphTitle = "Transcript-gene reconciliation"
		
	layout = dict(
				title = {'text': graphTitle, 'xanchor': 'auto', 'yanchor': 'auto'},
				  #paper_bgcolor='rgba(0,0,0,0)',
				  dragmode="select",
				  font=dict(family='Balto', size=10),
				  #width="100%",
				  height=800,
				  autosize=True,
				  showlegend=False,							   
				  xaxis=dict(
						range = [0, 100],
						titlefont=dict(
							family='Arial, sans-serif',
							size=15,
							color='black'
						),
						# animate=False,
						showticklabels=False,
						tickangle=45,
						tickfont=dict(
							family='Old Standard TT, serif',
							size=14,
							color='rgb(100,100,100)',
						),
						exponentformat='e',
						showexponent='all',					  
						showline=True,
						zeroline=False,
						showgrid=False,  # To visualize the vertical lines
						ticklen=4,
						#showticklabels=True,
						title=text_legend),
					yaxis=dict(
						# animate=False,
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
					paper_bgcolor='rgb(250,250,250)'
					#legend={'x': 0, 'y': 0}
				  )

	fig = go.Figure(data=nodes, layout=layout)
	return fig, [host_leaves_names, parasite_leaves_names]


def read_treefile(filename):
	tree = Phylo.read(filename, "newick")
	return tree


def read_tree_nw(tree_nw):
	tree = Phylo.read(StringIO(tree_nw), "newick")

	return tree
	
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

def textarea_example():
	file_tree = open("./Input/tmp_tree.xml", "r")
	tree_nw = file_tree.read()
	file_tree.close()
	return tree_nw


def get_recGeneSpecie_recProteinTree():
	tmp_file = open("./Input/tmp_tree.xml", 'rb')
	encoded_fig = base64.b64encode(tmp_file.read())
	tmp_file.close()
	decoded = base64.b64decode(encoded_fig)
	recTree = dataFromDoubleRecFile(StringIO(decoded.decode('utf-8')))
	if len(recTree) == 2 :
		if recTree[0][1] == "geneSpecie":
			recGeneSpecie = recTree[0][0]
			recProteinTree = recTree[1][0]
		elif recTree[0][1] == "transcriptGene":
			recGeneSpecie = recTree[1][0]
			recProteinTree = recTree[0][0]
	return recGeneSpecie, recProteinTree


def draw_example_figGeneSpecie():
	slanted = False
	recGeneSpecie, recProteinTree = get_recGeneSpecie_recProteinTree()
	figGeneSpecie, options_list2 = create_tree(recGeneSpecie, slanted, "geneSpecie", "red")
	return figGeneSpecie


def draw_example_figProteinGene():
	slanted = False
	recGeneSpecie, recProteinTree = get_recGeneSpecie_recProteinTree()
	figProteinGene, options_list = create_tree(recProteinTree, slanted, "transcriptGene", "blue")
	return figProteinGene


def get_local_fig_base64(path, format):
	tmp_file = open(path, 'rb')
	encoded_fig = base64.b64encode(tmp_file.read())
	tmp_file.close()
	return 'data:application/'+format+';base64,{}'.format(encoded_fig.decode())



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']			   
# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# webserver config
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, url_base_pathname='/DoubleRecViz/')
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

app.layout = html.Div(children=[
	html.H1(children='Double reconciliation visualization tool',style={'textAlign': 'center', 'color': '#000000'}),

	html.H6(children='''
		This web application allows visualizing transcript-gene-species tree reconciliations, as a combination of a gene tree embedded in a species tree and a transcript tree embedded inside a gene tree.
	''', style={'textAlign': 'center', 'color': '#000000'}),

	
	dcc.RadioItems(
		id = "inputdatamode",
		options=[
			{'label': 'Enter sequences data to compute tree', 'value': 'sequences'},
			{'label': 'Enter Reconciliated trees', 'value': 'tree'},
		],
		value='tree',
		style={'width': '100%','height':'50', 'textAlign': 'center'},
		labelStyle={'display': 'None'}
	),

	html.Div([		
		dcc.Dropdown(
			id='my-dropdown',
			options=[
				{'label': 'Slanted Tree', 'value': 'True'},
				{'label': 'No Slanted', 'value': 'False'},
			],
			style = {"display": "none"},
			value='slanted'
		),

		html.B("Gene sequences (fasta format):"),
		dcc.Textarea(
			id = "geneid",
			title = "",
			placeholder='>gene1 \nCGTATGGAATGCGTAAGAAGCAGGTCTGGGTAAGCATACGTGGGTAAGGGGAATGATTGAAAG\n>gene2\nCGTATGGAATGGGTAAGAAGCCAGTCTGGGTAAGCATACGTGGGTAAGGGGTTTGATTGAAAG',
			#value='This is a TextArea component',
			style={'width': '100%','height':'100'}
		),

		html.B("Transcript sequences (fasta format):"),
		dcc.Textarea(
			id = "transcriptid",
			title = "",
			placeholder='>transcript1\nATGCAAGCAGGTCTGGGGGAATGA\n>transcript2\nATGGAATGCAAGCAGCATACGTGGGGGAATGATTGA\n>transcript3\nATGGAAGCCAGTCTGGGGGTTTGA\n>transcript4\nATGGAATGGAAGCCACATACGTGGGGGTTTGATTGA',
			#value='>transcript1\nATGCAAGCAGGTCTGGGGGAATGA\n>transcript2\nATGGAATGCAAGCAGCATACGTGGGGGAATGATTGA\n>transcript3\nATGGAAGCCAGTCTGGGGGTTTGA\n>transcript4\nATGGAATGGAAGCCACATACGTGGGGGTTTGATTGA',
			style={'width': '100%','height':'100'}
		),

		html.B("Transcript to gene information:"),
		dcc.Textarea(
			id = "transcript2gene",
			title = "",
			placeholder='transcript1 gene1\ntranscript2 gene1\ntranscript3 gene2\ntranscript4 gene2',
			#value='transcript1 gene1\ntranscript2 gene1\ntranscript3 gene2\ntranscript4 gene2',
			style={'width': '100%','height':'100'}
		),

		html.B("Transcript exon position information:"),
		dcc.Textarea(
			id = "exonlist",
			title = "",
			placeholder='>transcript1\n0 4 8 12\n4 16 17 29\n16 24 48 56\n>transcript2\n0 9 3 12\n9 15 17 23\n15 24 34 43\n24 36 48 60\ntranscript3\n0 4 8 12\n4 16 17 29\n16 24 48 56\n>transcript4\n0 9 3 12\n9 15 17 23\n15 24 34 43\n24 36 48 60\n',
			#value='This is a TextArea component',
			style={'width': '100%','height':'100'}
		),
		html.B("Gene to species information:"),
		dcc.Textarea(
			id = "gene2species",
			title = "",
			placeholder='gene1 species1\ngene2 species2\ngene3 species3\ngene4 species3',
			style={'width': '100%','height':'100'}
		),		
		html.B("species tree:"),
		dcc.Textarea(
			id = "species tree",
			title = "",
			placeholder='(species2:1.00,((species1:0.45,species4:0.30):0.20,species3:0.32):0.22);',
			style={'width': '100%','height':'50'}
		),	
		html.Button('Computed Reconciliated trees', id='computed_trees', style={'textAlign': 'center', 'color': 'green', 'margin-left': '43%'}),
	],
	id ="sequences",
	style = {'display' : 'None'}
	),
	
	html.Div([
		dcc.Textarea(
			id = "proteinGeneSpecies", 
			title = "Protein Gene Reconciliation",   
			placeholder='Enter a value...',
			value=textarea_example(),
			style={'width': '100%','height':'150', 'margin-right': '9%'}
		),
		html.Button('Compute and draw', id='button', n_clicks=0, style={'textAlign': 'center', 'color': 'green', 'margin-left': '43%'}),
	],
	id ="trees",
	style = {'display' : 'None'}
	),
	dcc.ConfirmDialog(
			id='confirm',
			message='Input is not accepted, please check your input!!!',
		),
	dcc.ConfirmDialog(
			id='copute_error',
			message='Input Error: Please put <recGeneTree> before <recTransTree>',
		),
	dcc.Input(
		id='output-confirm',
		type='text',
		style={'display': 'none'}
		),
	html.Div(id= 'webserver_example_fig', children=[
		dcc.Graph(
		
			#style={
			 #   'width': 600,
			#	'height':600
			#},
			id='figGeneSpecie',
			#Gene Specie Reconciliation
			figure=draw_example_figGeneSpecie(),

			# Customize Download Plot figure's format:svg 
			config = {
			  'toImageButtonOptions': {
				'format': 'png', # one of png, svg, jpeg, webp
				'filename': 'figGeneSpecie',
				'height': 1000,
				'width': 2000,
				'scale': 1
			  }
			}
		),
		html.H6(children='Download Gene-species reconciliation result',style={'textAlign': 'left', 'color': '#000000'}),
		html.Div(children=[
			html.A(
			'figGeneSpecie.pdf',
			id='figGeneSpecie_pdf',
			download="figGeneSpecie.pdf",
			href=get_local_fig_base64("./Input/figGeneSpecie.pdf", "pdf"),
			target="_blank"
		)
		]),
		html.Div(children=[
		html.A(
			'figGeneSpecie.svg',
			id='figGeneSpecie_svg',
			download="figGeneSpecie.svg",
			href=get_local_fig_base64("./Input/figGeneSpecie.svg", "svg"),
			target="_blank"
		)
		]),
		dcc.Graph(
		
			#style={
			 #   'width': 600,
			#	'height':600
			#},
			id='figProteinGene',
			#Protein Gene Reconciliation',
			figure=draw_example_figProteinGene(),

			# Customize Download Plot figure's format:svg 
			config = {
			  'toImageButtonOptions': {
				'format': 'png', # one of png, svg, jpeg, webp
				'filename': 'figProteinGene',
				'height': 1000,
				'width': 2000,
				'scale': 1
			  }
			}
		),
		html.H6(children='Download Transcript-gene reconciliation result',style={'textAlign': 'left', 'color': '#000000'}),
		html.Div(children=[
		html.A(
			'figProteinGene.pdf',
			id='figProteinGene_pdf',
			download="figProteinGene.pdf",
			href=get_local_fig_base64("./Input/figProteinGene.pdf", "pdf"),
			target="_blank"
		)
		]),
		html.Div(children=[
		html.A(
			'figProteinGene.svg',
			id='figProteinGene_svg',
			download="figProteinGene.svg",
			href=get_local_fig_base64("./Input/figProteinGene.svg", "svg"),
			target="_blank"
		)
		]),
		]),
	html.Div(id='multiple_reconciliation_figGeneSpecie'),
])

"""
@app.callback(
	Output('proteinGene', 'value'),
	Output('geneSpecie', 'value'),
	[Input('computed_trees', 'n_clicks')]
)
def on_click(number_of_times_button_has_clicked):
	print("computed_trees")

@app.callback(
	Output('proteinGene', 'lang'),
	[Input('transcriptId', 'values'),
	Input('geneId', 'values'),
	Input('specieId', 'values')
	])
def set_cities_value(transcriptIdToRemove, GeneIdToRemove, specieIdToRemove):
	print(transcriptIdToRemove)
	print(GeneIdToRemove)
	print(specieIdToRemove)
"""       
@app.callback(
	Output('trees', 'style'),
	[Input('inputdatamode', 'value')])
def set_cities_value(inputdatamode):
	#print (inputdatamode)
	if inputdatamode == "tree":
		return {'display': 'block'}
	return dash.no_update


"""
@app.callback(
	Output('sequences', 'style'),
	[Input('inputdatamode', 'value')])
def set_cities_value(inputdatamode):
	#print (inputdatamode)
	if inputdatamode == "sequences":
		return {'display': 'block'}
	else:
		return {'display': 'None'}

""" 
"""			   
@app.callback(Output('figGeneSpecie', 'figure'),
	[Input('button', 'n_clicks'),
	Input('proteinGene', 'value'),
	Input('geneSpecie', 'value')])
def update_output(n_clicks, input2, input3):
	figGeneSpecie,options_list = create_tree(input3, False, "geneSpecie", "red")
	return figGeneSpecie
   

@app.callback(Output('figProteinGene', 'figure'),
	[Input('button', 'n_clicks'),
	Input('proteinGene', 'value'),
	Input('geneSpecie', 'value')])
def update_output(n_clicks, input2, input3):
	figProteinGene,options_list = create_tree(input2, False, "transcriptGene", "blue")
	return figProteinGene   
	
"""
@app.callback(Output('output-confirm', 'value'),
			  [Input('button', 'n_clicks')],
			  [State('proteinGeneSpecies', 'value')])
def display_confirm(n_clicks, input_xml):
	if n_clicks > 0:
		recPhylo_spTree_recGene = False
		recPhylo_gnTree_recTrans = False
		# tmp_xml = open("./Input/tmp_xml.xml", "w")
		# tmp_xml.write(input_xml)
		# tmp_xml.close()
		try:
			xml_file = lxml.etree.parse(StringIO(input_xml))
			xml_validator = lxml.etree.XMLSchema(file="./XSD/doubleRecPhylo.xsd")
			if xml_file.getroot().tag == "recPhylo":
				if xml_file.find(".//spTree"):
					recPhylo_spTree_recGene = True
					xml_validator = lxml.etree.XMLSchema(file="./XSD/recPhylo_spTree_recGene.xsd")
				if xml_file.find(".//gnTree"):
					recPhylo_gnTree_recTrans = True
					xml_validator = lxml.etree.XMLSchema(file="./XSD/recPhylo_gnTree_recTrans.xsd")
			is_valid = xml_validator.validate(xml_file)
			if is_valid:
				if recPhylo_spTree_recGene:
					return "YES_recPhylo_spTree_recGene"
				if recPhylo_gnTree_recTrans:
					return "YES_recPhylo_gnTree_recTrans"
				return "YES"
			return "NO"
		except lxml.etree.XMLSyntaxError:
			return "NO"
	return dash.no_update


@app.callback(Output('confirm', 'displayed'),
			  [Input('output-confirm', 'value'), Input('button', 'n_clicks')])
def display_confirm(value, n_clicks):
	if n_clicks > 0:
		if value == 'NO':
			return True
		return False
	return dash.no_update



@app.callback(
	[
	Output('multiple_reconciliation_figGeneSpecie','children'),
	Output('webserver_example_fig', 'style'),
	Output('copute_error', 'displayed'),
	],
	[Input('button', 'n_clicks'), Input('output-confirm', 'value')],
	[State('proteinGeneSpecies', 'value')],
	)
def update_output(n_clicks, output_confirm, input2):
	if n_clicks > 0:
		display_error_message = False
		if output_confirm in ["YES", "YES_recPhylo_spTree_recGene", "YES_recPhylo_gnTree_recTrans"]:
			# add multiple reconciliation figs
			figGeneSpecie_figProteinGene = dict()
			error_message = ""

			multiple_reconciliation_figGeneSpecie_list = []
			multiple_reconciliation_figs = []

			if output_confirm == "YES":
				# tmp_tree = open("./Input/tmp_xml.xml", "w")
				# tmp_tree.write(input2)
				# tmp_tree.close()
				recTree = dataFromDoubleRecFile(StringIO(input2))
				
				for i in range(len(recTree)):
					if recTree[i][1] == "geneSpecie":
						fig_GeneSpecie,options_list = create_tree(recTree[i][0], False, "geneSpecie", "red")
						multiple_reconciliation_figGeneSpecie_list.append(fig_GeneSpecie)

					if recTree[i][1] == "transcriptGene":
						fig_ProteinGene,options_list = create_tree(recTree[i][0], False, "transcriptGene", "blue")
						if multiple_reconciliation_figGeneSpecie_list:
							key = len(multiple_reconciliation_figGeneSpecie_list) - 1
							if key not in figGeneSpecie_figProteinGene:
								figGeneSpecie_figProteinGene[key] = []
							figGeneSpecie_figProteinGene[key].append(fig_ProteinGene)
						else:
							display_error_message = True
						
				for i in range(len(multiple_reconciliation_figGeneSpecie_list)):
					iteration_figGeneSpecie = multiple_reconciliation_figGeneSpecie_list[i]
					iteration_figProteinGenes = figGeneSpecie_figProteinGene[i]

					# encoded_figGeneSpecie_pdf = base64.b64encode(iteration_figGeneSpecie.to_image(format="pdf", engine="kaleido", width=2000, height=1000, scale=2))
					# encoded_figGeneSpecie_svg = base64.b64encode(iteration_figGeneSpecie.to_image(format="svg", engine="kaleido", width=2000, height=1000, scale=2))
					encoded_figGeneSpecie_pdf = base64.b64encode(iteration_figGeneSpecie.to_image(format="pdf", width=2000, height=1000, scale=2))
					encoded_figGeneSpecie_svg = base64.b64encode(iteration_figGeneSpecie.to_image(format="svg", width=2000, height=1000, scale=2))
					div = html.Div([
						dcc.Graph(figure=iteration_figGeneSpecie, config = {'toImageButtonOptions': {'format': 'png', 'filename': 'figGeneSpecie', 'height': 1000, 'width': 2000, 'scale': 1}}),
						html.H6(children='Download Gene-species reconciliation result',style={'textAlign': 'left', 'color': '#000000'}),
						html.Div(children=[
							html.A(
							'figGeneSpecie.pdf',
							download="figGeneSpecie.pdf",
							href='data:application/pdf;base64,{}'.format(encoded_figGeneSpecie_pdf.decode()),
							target="_blank"
							)
						]),
						html.Div(children=[
						html.A(
							'figGeneSpecie.svg',
							download="figGeneSpecie.svg",
							href='data:application/svg;base64,{}'.format(encoded_figGeneSpecie_svg.decode()),
							target="_blank"
							)
						])
					])
					multiple_reconciliation_figs.append(div)
					for recProteinTree in iteration_figProteinGenes:
						# encoded_figGeneSpecie_pdf = base64.b64encode(recProteinTree.to_image(format="pdf", engine="kaleido", width=2000, height=1000, scale=2))
						# encoded_figGeneSpecie_svg = base64.b64encode(recProteinTree.to_image(format="svg", engine="kaleido", width=2000, height=1000, scale=2))
						encoded_figGeneSpecie_pdf = base64.b64encode(recProteinTree.to_image(format="pdf", width=2000, height=1000, scale=2))
						encoded_figGeneSpecie_svg = base64.b64encode(recProteinTree.to_image(format="svg", width=2000, height=1000, scale=2))
						div = html.Div([
							dcc.Graph(figure=recProteinTree, config = {'toImageButtonOptions': {'format': 'png', 'filename': 'figGeneSpecie', 'height': 1000, 'width': 2000, 'scale': 1}}),
							
							html.H6(children='Download Transcript-gene reconciliation result',style={'textAlign': 'left', 'color': '#000000'}),
							html.Div(children=[
								html.A(
								'figProteinGene.pdf',
								download="figProteinGene.pdf",
								href='data:application/pdf;base64,{}'.format(encoded_figGeneSpecie_pdf.decode()),
								target="_blank"
							)
							]),
							html.Div(children=[
							html.A(
								'figProteinGene.svg',
								download="figProteinGene.svg",
								href='data:application/svg;base64,{}'.format(encoded_figGeneSpecie_svg.decode()),
								target="_blank"
							)
							]),
						])
						multiple_reconciliation_figs.append(div)


				
			if output_confirm == "YES_recPhylo_spTree_recGene":
				first_figGeneSpecie,options_list = create_tree(input2, False, "geneSpecie", "red")
				# encoded_figGeneSpecie_pdf = base64.b64encode(first_figGeneSpecie.to_image(format="pdf", engine="kaleido", width=2000, height=1000, scale=2))
				# encoded_figGeneSpecie_svg = base64.b64encode(first_figGeneSpecie.to_image(format="svg", engine="kaleido", width=2000, height=1000, scale=2))
				encoded_figGeneSpecie_pdf = base64.b64encode(first_figGeneSpecie.to_image(format="pdf", width=2000, height=1000, scale=2))
				encoded_figGeneSpecie_svg = base64.b64encode(first_figGeneSpecie.to_image(format="svg", width=2000, height=1000, scale=2))
				div = html.Div([
					dcc.Graph(figure=first_figGeneSpecie, config = {'toImageButtonOptions': {'format': 'png', 'filename': 'figGeneSpecie', 'height': 1000, 'width': 2000, 'scale': 1}}),
					html.H6(children='Download Gene-species reconciliation result',style={'textAlign': 'left', 'color': '#000000'}),
					html.Div(children=[
						html.A(
						'figGeneSpecie.pdf',
						download="figGeneSpecie.pdf",
						href='data:application/pdf;base64,{}'.format(encoded_figGeneSpecie_pdf.decode()),
						target="_blank"
						)
					]),
					html.Div(children=[
					html.A(
						'figGeneSpecie.svg',
						download="figGeneSpecie.svg",
						href='data:application/svg;base64,{}'.format(encoded_figGeneSpecie_svg.decode()),
						target="_blank"
						)
					])
				])
				multiple_reconciliation_figs.append(div)

			if output_confirm == "YES_recPhylo_gnTree_recTrans":
				first_recProteinTree,options_list = create_tree(input2, False, "transcriptGene", "blue")
				# encoded_figProteinGene_pdf = base64.b64encode(first_recProteinTree.to_image(format="pdf", engine="kaleido", width=2000, height=1000, scale=2))
				# encoded_figProteinGene_svg = base64.b64encode(first_recProteinTree.to_image(format="svg", engine="kaleido", width=2000, height=1000, scale=2))
				encoded_figProteinGene_pdf = base64.b64encode(first_recProteinTree.to_image(format="pdf", width=2000, height=1000, scale=2))
				encoded_figProteinGene_svg = base64.b64encode(first_recProteinTree.to_image(format="svg", width=2000, height=1000, scale=2))
				div = html.Div([
					dcc.Graph(figure=first_recProteinTree, config = {'toImageButtonOptions': {'format': 'png', 'filename': 'figGeneSpecie', 'height': 1000, 'width': 2000, 'scale': 1}}),
					
					html.H6(children='Download Transcript-gene reconciliation result',style={'textAlign': 'left', 'color': '#000000'}),
					html.Div(children=[
						html.A(
						'figProteinGene.pdf',
						download="figProteinGene.pdf",
						href='data:application/pdf;base64,{}'.format(encoded_figGeneSpecie_pdf.decode()),
						target="_blank"
					)
					]),
					html.Div(children=[
					html.A(
						'figProteinGene.svg',
						download="figProteinGene.svg",
						href='data:application/svg;base64,{}'.format(encoded_figGeneSpecie_svg.decode()),
						target="_blank"
					)
					]),
				])
				multiple_reconciliation_figs.append(div)

			result_multiple_reconciliation = html.Div(children = multiple_reconciliation_figs)

			return result_multiple_reconciliation, {'display' : 'None'}, display_error_message
	return dash.no_update


import argparse
def build_arg_parser():
	parser = argparse.ArgumentParser(description="Run DoubleRecViz webserver")
	parser.add_argument('-d', '--debug', help="debug output option: True or False default(True)", default="True")
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	args = parser.parse_args()
	d = args.debug
	display_debug = True
	if(d == "False"):
		display_debug = False
	app.run_server(debug=display_debug)
