import networkx as nx
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from scipy.spatial import distance
from itertools import count

#1. Read Data
D = []
with open("norm_cellcycle_384_17.txt") as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        D.append(line[:19])

D = D[1:]
n = len(D)
#n=5
#2. Create a graph
G = nx.Graph()

#2.1 Define nodes and attributes 
G.add_nodes_from([(i,{"label": D[i][0],"group":int(D[i][1]), "data":D[i][2:]} ) for i in range(n)])
print(len(G.nodes()))

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

#2.3 Add degree 


#nx.set_node_attributes(G, Degree, name='degree')


#nx.draw(G,with_labels=True,  font_size=5, font_color= 'white',  edge_color="grey", node_size=[v * 10 for v in Degree.values()], node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#nx.draw(G,  with_labels=True,  font_size=5, font_color= 'white',edge_color="grey", node_color=[G.nodes[i]["group"] for i in range(n)], pos=nx.fruchterman_reingold_layout(G))
#plt.show()
#3. Draw Graph accdg to color