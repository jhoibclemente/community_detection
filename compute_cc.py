import networkx as nx
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from scipy.spatial import distance
from itertools import count
import hvplot.networkx as hvnx
from networkx.algorithms import approximation


delta = float(sys.argv[1])
print(delta)
\\3D = []
with open("norm_cellcycle_384_17.txt") as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        D.append(line[:19])

D = D[1:]
n = len(D)
#n = 20
#2. Create a graph
G = nx.Graph()

#2.1 Define nodes and attributes 
G.add_nodes_from([(i,{"label": D[i][0],"group":int(D[i][1]), "data":D[i][2:]} ) for i in range(n)])
print("Number of Nodes", len(G.nodes()))
print("Number of Edges", G.number_of_edges())
print("Transitivity", nx.transitivity(G))


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
#3. Compute for CC
CC = [G.subgraph(c).copy() for c in nx.connected_components(G) if len(c)>1]
for c  in CC:
    print(len(c), list(c.nodes(data="label")))
#print([[len(c), [node.nodes for node in c]] for c in sorted(nx.connected_components(G), key=len, reverse=True) if len(c)>1])
#nx.draw(G,with_labels=True,  font_size=5, font_color= 'white',  edge_color="grey", node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#nx.draw(G, with_labels=True,  font_size=5, font_color='white',edge_color="grey", node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#plt.show()
#img = hvnx.draw(G, with_labels=False, node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#hvnx.show(img)
