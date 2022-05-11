import networkx as nx
from sklearn.metrics import rand_score
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import f1_score
import csv
import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from scipy.spatial import distance
from itertools import count
import hvplot.networkx as hvnx
from networkx.algorithms import approximation
from networkx.algorithms import community
import community.community_louvain as community_louvain
import plotly.graph_objects as go
import pandas as pd

def spring_3d_viz():
    spring_3D = nx.fruchterman_reingold_layout(G, dim =3)
    #we need to seperate the X,Y,Z coordinates for Plotly
    x_nodes = [spring_3D[i][0] for i in range(n)]# x-coordinates of nodes
    y_nodes = [spring_3D[i][1] for i in range(n)]# y-coordinates
    z_nodes = [spring_3D[i][2] for i in range(n)]# z-coordinates

    #We also need a list of edges to include in the plot
    edge_list = G.edges()

    edge_list
    #we  need to create lists that contain the starting and ending coordinates of each edge.
    x_edges=[]
    y_edges=[]
    z_edges=[]

    #need to fill these with all of the coordiates
    for edge in edge_list:
        #format: [beginning,ending,None]
        x_coords = [spring_3D[edge[0]][0],spring_3D[edge[1]][0],None]
        x_edges += x_coords

        y_coords = [spring_3D[edge[0]][1],spring_3D[edge[1]][1],None]
        y_edges += y_coords

        z_coords = [spring_3D[edge[0]][2],spring_3D[edge[1]][2],None]
        z_edges += z_coords
    #create a trace for the edges
    trace_edges = go.Scatter3d(x=x_edges,
                            y=y_edges,
                            z=z_edges,
                            mode='lines',
                            line=dict(color='grey', width=0.5),
                            hoverinfo='none')
    #create a trace for the nodes
    trace_nodes = go.Scatter3d(x=x_nodes,
                            y=y_nodes,
                            z=z_nodes,
                            mode='markers',
                            marker=dict(symbol='circle',
                                        size=[v * 4 for v in Degree.values()],
                                        color=list(nx.get_node_attributes(G,'group').values()), #color the nodes according to their community
                                        colorscale=['purple', 'green', 'blue', 'pink', 'yellow'], #either green or mageneta
                                        line=dict(color='black', width=0.1)),
                                        text =list(nx.get_node_attributes(G,'label').values()))
                                        
                            #text=club_labels,
                            #hoverinfo='label'
    #we need to set the axis for the plot 
    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title='')                       

    #also need to create the layout for our plot
    layout = go.Layout(title="Community Detection of Yeast's Gene Expression Data using delta = "+str(delta),
                    width=1200,
                    height=800,
                    showlegend=False,
                    scene=dict(xaxis=dict(axis),
                            yaxis=dict(axis),
                            zaxis=dict(axis),
                            ),
                    margin=dict(t=100),
                    hovermode='closest')
    
    data = [trace_edges, trace_nodes]
    fig = go.Figure(data=data, layout=layout)

    fig.show()

delta = float(sys.argv[1])
print(delta)
D = []
with open("norm_cellcycle_384_17.txt") as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        D.append(line[:19])

D = D[1:]
n = len(D)
n = 20
#2. Create a graph
G = nx.Graph()

#2.1 Define nodes and attributes 
G.add_nodes_from([(i,{"label": D[i][0],"group":int(D[i][1]), "data":D[i][2:]} ) for i in range(n)])



#2.2 Define edges 
for i in range(n):
    for j in range(i+1, n):
  
        #print(G.nodes[i])
        #print(G.nodes[j])
        x = [float(xi) for xi in G.nodes[i]["data"]]
        y = [float(yi) for yi in G.nodes[j]["data"]]
     
        d = distance.euclidean(x,y)
        m = distance.minkowski(x,y)
        pearson = np.corrcoef(x,y)[0][1]
        spearman, ps = spearmanr(x,y)
        kendall, pk = kendalltau(x,y)
        if pearson>=delta:
            G.add_edge(i,j, euc = d, mink= m, weight = pearson, spearman = spearman, kendall = kendall) 
print(approximation.average_clustering(G, trials=1000, seed=10))
#2.3 Add node attributes (print properties)

Degree = dict(G.degree)
nx.set_node_attributes(G, Degree, 'degree')
print("Number of Nodes", len(G.nodes()))
print("Number of Edges", G.number_of_edges())
print("Clustering Coeff:", approximation.average_clustering(G))

#3. Compute for CC
CC = [G.subgraph(c).copy() for c in nx.connected_components(G) if len(c)>1]
max = -1
BigC = []
for c  in CC:
    print(len(c), list(c.nodes(data="label")))
    if len(c) > max: 
        max = len(c)
        BigC = c 
#4.1 Cluster Biggest  Connected Component


print("Max Component with", len(BigC))
print("Clustering Coeff: ", approximation.average_clustering(BigC))
nx.set_node_attributes(G, -1, 'gn_grp')

Greedy = community.greedy_modularity_communities(BigC, n_communities=5)

x = [list(sorted(i)) for i in Greedy]
print(x)


gn_dic = {}
i = 0
for c in x:
    print("Group", i, c)
    for j in c: 
        #print(j, i)
        gn_dic[j] = i
    i = i+1

nx.set_node_attributes(G, gn_dic, 'gn_grp')
#4.2 Adjusted Rand Index Computation
A = [G.nodes[i]["group"] for i in range(n)]
B = [G.nodes[i]["gn_grp"] for i in range(n)]
GroundT = []
Predicted = []
for i in range(n):
    if B[i] != -1:
        GroundT.append(A[i])
        Predicted.append(B[i])

print(len(GroundT)==len(Predicted))
print("Rand Score:", rand_score(Predicted, GroundT))
print("Adjusted Rand Score:", adjusted_rand_score(Predicted, GroundT))
print("F1 Score:", f1_score(Predicted, GroundT, average='micro'))

print("True: ", GroundT)
print("Predicted: ", Predicted) 


#5. Visualization
#option 1
#nx.draw(G,with_labels=True,  font_size=5, font_color= 'white',  edge_color="grey",  node_color=[G.nodes[i]["gn_grp"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#plt.show()

#option 2
#nx.draw(G,with_labels=True,  font_size=5, font_color= 'white',  edge_color="grey", node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i]["gn_grp"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#plt.show()

#option 3
#img = hvnx.draw(G, with_labels=False, node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#hvnx.show(img)

#option 4 (removed isolated nodes)
#img = hvnx.draw(G, edge_color="grey", with_labels=False, node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#hvnx.show(img)

#option 5 3d
#spring_3d_viz()