# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from Bio import Phylo
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import random
from dash.dependencies import Input, Output

app = dash.Dash()

def label_tree_internal_nodes(tree):
    i = 1
    for t in tree.find_clades():
        if t.name == None:
            t.name = str(i)
            i = i + 1

def get_x_coordinates_new(tree, width):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
        
    for clade in xcoords.keys():
        if len(clade.clades) ==2:
            xcoords[clade] = [xcoords[clade], xcoords[clade], xcoords[clade], xcoords[clade], xcoords[clade]+width, xcoords[clade]+ width]
        else:
            xcoords[clade] = [xcoords[clade], xcoords[clade]]
			
    return xcoords
    
    
def get_x_coordinates(tree):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords

def get_y_coordinates_new(tree, width, dist=1.3):
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
            y_child_bottom = y_coords[child_bottom]
            y_child_top = y_coords[child_top]        	        	
            ycoords[clade] = [y_child_bottom, ycoords[clade], ycoords[clade] + width, y_child_top, y_child_bottom + width,  y_child_top - width]        
        else: 
            ycoords[clade] = [ycoords[clade], ycoords[clade] + width]        
            
    return ycoords
    
def get_y_coordinates(tree, dist=1.3):
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
    return ycoords

def get_clade_lines_new(line_shapes, width, orientation='horizontal',y_curr2=[0,0], y_curr=[0,0], x_start=[0,0], x_curr=[0,0], y_bot=[0,0], y_top=[0,0], line_color='rgb(25,25,25)', line_width=3):
    """define a shape of type 'line', for branch
    """

    branch_line_1 = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    branch_line_2 = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )             
    branch_line_3 = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )             
    branch_line_4 = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )                                                                     
    if orientation == 'horizontal':
    
    	if y_curr[0] > 
        branch_line_1.update(x0=x_start[0] + width,
                           y0=y_curr[0],
                           x1=x_curr[0],
                           y1=y_curr[0])
        line_shapes.append(branch_line_1) 
                          
        branch_line_2.update(x0=x_start[1],
                           y0=y_curr[1],
                           x1=x_curr[1],
                           y1=y_curr[1])
        line_shapes.append(branch_line_2) 
                           
        print (orientation," x0 = ",x_start,  " y0 = ",y_curr, " x1 = ", x_curr, " y1 = ",y_curr)
    elif orientation == 'vertical':
        branch_line_1.update(x0=x_curr[0],
                           y0=y_curr2[0],
                           x1=x_curr[0],
                           y1=y_bot[0])
        line_shapes.append(branch_line_1) 
                      
        branch_line_2.update(x0=x_curr[1],
                           y0=y_curr2	[1],
                           x1=x_curr[1],
                           y1=y_top[1])                           
        line_shapes.append(branch_line_2) 
		
        branch_line_3.update(x0=x_curr[0] + width,
                           y0=y_top[0],
                           x1=x_curr[0] + width,
                           y1=y_bot[1])
        line_shapes.append(branch_line_3) 
                
        print (orientation, " x_curr = ",x_curr,  " y_bot = ",y_bot, " x_curr = ", x_curr, " y_top = ",y_top)                           
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")
   
    print ("------------------------------")
    return line_shapes


def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=3):
    """define a shape of type 'line', for branch
    """
    #print(y_curr, x_start, x_curr, y_bot, y_top)

    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
        print (orientation," x0 = ",x_start,  " y0 = ",y_curr, " x1 = ", x_curr, " y1 = ",y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
        print (orientation, " x0 = ",x_curr,  " y0 = ",y_bot, " x1 = ", x_curr, " y1 = ",y_top)                           
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")
   
    print ("------------------------------")
    return branch_line


def draw_clade_new(clade, x_start, line_shapes, width, line_color='rgb(15,15,15)', line_width=3, x_coords=[0,0], y_coords=[0,0]):
    """Recursively draw the tree branches, down from the given clade"""    
    print("nombre fils", len(clade.clades))
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]
    print(x_curr, y_curr)
    # Draw a horizontal line from start to here
    line_shapes = get_clade_lines_new(line_shapes, width, orientation='horizontal', y_curr2=y_curr, y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)

    #line_shapes.append(branch_line)
   
    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[-1]]
        y_bot = y_coords[clade.clades[0]]
        print(len(clade.clades), y_top, y_bot)
        line_shapes = (get_clade_lines_new(line_shapes, width, orientation='vertical',y_curr2=y_curr, x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade_new(child, x_curr, line_shapes, width, posistion= x_coords=x_coords, y_coords=y_coords)
            

def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=3, x_coords=0, y_coords=0):
    """Recursively draw the tree branches, down from the given clade"""    
    print("nombre fils", len(clade.clades))
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]
    print(x_curr, y_curr)
    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)

    line_shapes.append(branch_line)
   
    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]
        print(len(clade.clades), y_top, y_bot)
        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)

            
def random_color():
    # Make a list of colors to picvk from
    colors = ["red", "green", "blue", "orange", "purple", "pink", "yellow", "black", "gray", "sliver", "violet",
              "yellowgreen", "turquoise", "sienna", "salmon"]
    color = random.choice(colors)
    return tuple(color)
                
def create_tree(tree_file):
    width= 0.2
    tree = read_treefile(tree_file) 
    label_tree_internal_nodes(tree)  
    x_coords_new = get_x_coordinates_new(tree)
    y_coords_new = get_y_coordinates_new(tree, width, 1.3)
        
    for clade in x_coords_new.keys():
         print(clade, x_coords_new[clade])
         
    for clade in y_coords_new.keys():
         print(clade, y_coords_new[clade])         
 
    #x_coords = get_x_coordinates(tree)
    #y_coords = get_y_coordinates(tree)
    #print(x_coords)
    #print(y_coords)
    
    line_shapes = []
    #draw_clade(tree.root, 0, line_shapes, line_color=random_color(), line_width=3, x_coords=x_coords, y_coords=y_coords)
    draw_clade_new(tree.root, [0,0], line_shapes, 20*width, line_color=random_color(), line_width=3, x_coords=x_coords_new,
               y_coords=y_coords_new)

               
    my_tree_clades = x_coords_new.keys()
    X = []
    Y = []
    text = []

    for cl in my_tree_clades:
        X.append(x_coords_new[cl][0])
        Y.append(y_coords_new[cl][0])
        X.append(x_coords_new[cl][1])
        Y.append(y_coords_new[cl][1])        
        text.append(cl.name)
    print("=========================================================")
    print(X, Y, text)
    print(line_shapes)
    print("=========================================================")    

    """
    df = read_metadata(metadata_file)
    data_metadata_stat_csv = df.groupby(ord_by)['Strain'].count()

    # for index_val, series_val in data_metadata_stat_csv.iteritems():
    df.columns
    nb_genome = len(df)

    graph_title = create_title(virus_name, nb_genome)
    """
    intermediate_node_color = 'rgb(100,100,100)'


    country = []
    region = []
    color = [intermediate_node_color] * len(X)

    

    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
                )

    label_legend = ["a", "b", "c"]
    nodes = []
    
    for elt in label_legend:
        node = dict(type='scatter',
                x=X,
                y=Y,
                mode='markers',
                marker=dict(color=color, size=5),
                text=text,  # vignet information of each node
                hoverinfo='hoverinfo',
                name=elt
                )
        nodes.append(node)
    print(nodes)
    layout = dict(title="Titre Esaie",
                  paper_bgcolor='rgba(0,0,0,0)',
                  dragmode="select",
                  font=dict(family='Balto', size=14),
                  width=1000,
                  height=1000,
                  #autosize=True,
                  showlegend=True,
                  xaxis=dict(showline=True,
                             zeroline=False,
                             showgrid=True,  # To visualize the vertical lines
                             ticklen=4,
                             showticklabels=True,
                             title='branch length'),
                  yaxis=axis,
                  hovermode='closest',
                  shapes=line_shapes,
                  plot_bgcolor='rgb(250,250,250)',
                  legend={'x': 0, 'y': 1},
                  )

    fig = dict(data=nodes, layout=layout)
    print(layout)
    return fig


def read_treefile(filename):
    tree = Phylo.read(filename, "newick")
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

    


fig = create_tree("/home/local/USHERBROOKE/kuie2201/Documents/Dash/TreeVisualization/test.nw")
                

app.layout = html.Div(children=[
    html.H1(children='Reconciliation visualisation tool'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dcc.Graph(
        id='graph-1',
        figure=fig
    )
])


if __name__ == '__main__':
    app.run_server(debug=True)
