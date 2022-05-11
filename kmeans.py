import networkx as nx
from sklearn.metrics import rand_score
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import f1_score
from scipy.cluster.hierarchy import dendrogram
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
from cdlib import algorithms as alg
from networkx.algorithms import approximation
from networkx.algorithms import community
import community.community_louvain as community_louvain
import plotly.graph_objects as go
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN
from sklearn.cluster import Birch
from sklearn.cluster import AgglomerativeClustering

k = int(sys.argv[1])
D= []
with open("input/TEMP_PH.csv") as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        D.append(line[:19])

D = D[1:]
n = len(D)
#n = 20

true_labels = [(D[i][0]) for i in range(n)]
features = np.array([D[i][5:] for i in range(n)])


clustering = KMeans(n_clusters=k, init='k-means++',random_state=0, algorithm='elkan').fit(features)
#clustering = SpectralClustering(n_clusters=k,random_state=0, assign_labels='discretize').fit(features)
#clustering = DBSCAN(eps=3,min_samples=k).fit(features)

#clustering = Birch(n_clusters=k).fit(features)
#Predicted = list(clustering.predict(features))

#clustering = AgglomerativeClustering(n_clusters=k).fit(features)

Predicted = list(clustering.labels_)



GroundT = true_labels
print(Predicted)
print(GroundT)

#print("k, Rand, ARI, F1")
#print(rand_score(Predicted, GroundT), ', ',adjusted_rand_score(Predicted, GroundT), ', ',f1_score(Predicted, GroundT, average='micro'))
