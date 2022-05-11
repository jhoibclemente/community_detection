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
import plotly.graph_objects as go
import pandas as pd

def spring_3d_viz(G, coloring):
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
                                        #size=[v * 1 for v in Degree.values()],
                                        size = 5,
                                        color=list(nx.get_node_attributes(G,coloring).values()), #color the nodes according to their community
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

delta = int(sys.argv[1])
#k = int(sys.argv[2])
#coloring = str(sys.argv[2])
data = str(sys.argv[2])
coloring = "group"
print(delta)
D = []
with open(data) as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        D.append(line[:19])

D = D[1:]
n = len(D)
#n = 20
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
        G.add_edge(i,j, euc = d, mink= m, weight = pearson, spearman = spearman, kendall = kendall) 
        #G.add_edge(i,j, weight=spearman)
        #print(i,j, spearman)

print(approximation.average_clustering(G, trials=1000, seed=10))

#Remove Extra Edges
NG= nx.Graph()
NG.add_nodes_from([(i,{"label": D[i][0],"group":int(D[i][1]), "data":D[i][2:]} ) for i in range(n)])
#NG.add_nodes_from(G.nodes())
for i in range(n):
    neighbors = list(sorted(G[i].items(), reverse=True,key=lambda edge: edge[1]['weight']))[:delta]
    #print(i,[(i, node) for node, _metadata in neighbors])
    
    for (j,data) in neighbors: NG.add_edge(i,j, weight=data['weight'])
 
#2.3 Add node attributes (print properties)

Degree = dict(NG.degree)
nx.set_node_attributes(G, Degree, 'degree')
#print("Number of Nodes", len(G.nodes()))
#print("Number of Edges", G.number_of_edges())
#print("Clustering Coeff:", approximation.average_clustering(G))

#3. Compute for CC
CC = [NG.subgraph(c).copy() for c in nx.connected_components(NG) if len(c)>1]
max = -1
BigC = []
for c  in CC:
    #print(len(c), list(c.nodes(data="group")))
    if len(c) > max: 
        max = len(c)
        BigC = c 
#4.1 Cluster Biggest  Connected Component


print("Max Component with", len(BigC))
print("Clustering Coeff: ", approximation.average_clustering(BigC))
nx.set_node_attributes(NG, -1, 'gn_grp')
#Option 1 Girvan Newman
for k in range(1,11):

    CGN = community.girvan_newman(BigC)
    limited = itertools.takewhile(lambda c: len(c) <= k, CGN) 
    cgn = nx.Graph()
    for c in limited:
        if len(c)== k: cgn = c

    x = [list(sorted(i)) for i in cgn]
    #print(x)


    gn_dic = {}
    i = 0
    for c in x:
        #print("Group", i, c)
        for j in c: 
            #print(j, i)
            gn_dic[j] = i
        i = i+1

    nx.set_node_attributes(NG, gn_dic, 'gn_grp')
    #4.2 Adjusted Rand Index Computation
    A = [NG.nodes[i]["group"] for i in range(n)]
    B = [NG.nodes[i]["gn_grp"] for i in range(n)]
    GroundT = []
    Predicted = []
    for i in range(n):
        if B[i] != -1:
            GroundT.append(A[i])
            Predicted.append(B[i])

    
    #print(len(GroundT)==len(Predicted))
    #print(k, "Rand, ARI, F1")
    print(delta,', ', k,', ', rand_score(Predicted, GroundT), ', ',adjusted_rand_score(Predicted, GroundT), ', ',f1_score(Predicted, GroundT, average='micro'))


#5. Visualization
#option 1
#nx.draw(G,with_labels=True,  font_size=5, font_color= 'white',  edge_color="grey",  node_color=[G.nodes[i]["gn_grp"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#plt.show()

#option 2
#nx.draw(G,with_labels=True,  font_size=5, font_color= 'white',  edge_color="grey", node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i]["gn_grp"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#plt.show()

#option 3
#img1 = hvnx.draw(G, title="Rank based d="+str(delta),with_labels=False, node_color=[G.nodes[i][coloring] for i in range(n)], node_size=[v * 10 for v in Degree.values()],  pos=nx.fruchterman_reingold_layout(G))
#hvnx.show(img1)

#option 4 (removed isolated nodes)
#img = hvnx.draw(G, edge_color="grey", with_labels=False, node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i][coloring] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#hvnx.show(img)

#option 5 3d
#spring_3d_viz(G,"group")

#option 6 planar non trivial cc
#img = hvnx.draw(G, edge_color="grey", with_labels=False, node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i][coloring] for i in range(n)], pos=nx.planar_layout(G))
#hvnx.show(img)