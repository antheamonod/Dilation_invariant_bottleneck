# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 08:50:07 2021

@author: Administrator
"""

import numpy as np
import matplotlib.pyplot as plt
import DI_distance as dist
from sklearn import manifold
import umap

# load data
    
file_list = ['activity_net_euclidean_cosine.npz','activity_net_euclidean.npz',
             'activity_net_euclidean_expmap_poincare.npz','activity_net_poincare.npz',
             'mammals_euclidean_cosine.npz','mammals_euclidean.npz',
             'mammals_euclidean_expmap_poincare.npz','mammals_poincare.npz']


dist_list = [np.load(u)['arr_0'] for u in file_list]
Dmatrix = np.zeros((8,8))
 
for i in range(0,8):
    for j in range(0,8):
        dataA = np.delete(dist_list[i][0], -1, axis = 1)
        dataB = np.delete(dist_list[j][0], -1, axis = 1)
        ep = 1e-10
        dataA = np.log(dataA+ep)
        dataB = np.log(dataB+ep)
        dilaInv, bottleneckDistance, minBottleneckDistance = dist.myDistance(dataA, dataB, 10)
        Dmatrix[i,j] = minBottleneckDistance

        
# mds visualization
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=10,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(Dmatrix).embedding_


plt.figure()
plt.scatter(pos[0:2, 0], pos[0:2, 1], color='red', label='Ac-Eu')
plt.scatter(pos[2:4, 0], pos[2:4, 1], color='blue', label='Ac-Po')
plt.scatter(pos[4:6, 0], pos[4:6, 1], color='green', label='Ma-Eu')
plt.scatter(pos[6:8, 0], pos[6:8, 1], color='black', label='Ma-Po')
label = ['acNet-euclidean','acNet-poincare','mammal-euclidean','mammal-poincare']
plt.legend(label)


# umap visualization
RANDOM_SEED = 1 
reducer = umap.UMAP(
    n_neighbors = 4,       # default: 15
    n_components = 2,       # 2D atlas
    metric = 'precomputed', # we have already computed the pairwise distances
    min_dist = .05,         # default: 0.1
    spread = 1,             # default: 1
    random_state = RANDOM_SEED)

embedding = reducer.fit_transform(Dmatrix) 

# Plot the UMAP atlas 
plt.figure()
plt.scatter(embedding[0:2, 0], embedding[0:2, 1], color='red', label='Ac-Eu')
plt.scatter(embedding[2:4, 0], embedding[2:4, 1], color='blue', label='Ac-Po')
plt.scatter(embedding[4:6, 0], embedding[4:6, 1], color='green', label='Ma-Eu')
plt.scatter(embedding[6:8, 0], embedding[6:8, 1], color='black', label='Ma-Po')
label = ['acNet-euclidean','acNet-poincare','mammal-euclidean','mammal-poincare']
plt.legend(label)

