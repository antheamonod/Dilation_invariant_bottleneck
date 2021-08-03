# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 09:16:20 2021

@author: Administrator
"""

import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import DI_distance as dist
from sklearn import manifold
import umap
from gudhi.wasserstein.barycenter import lagrangian_barycenter

def read_pickle_object(path):
    with open(path, 'rb') as handle:
        b = pickle.load(handle)
    return b

def make_diagramlist(pathname, diagramlist):
    path = os.getcwd()+"\\PH Diagrams\\" + pathname.replace(' ','_') + "\\"
    pklfiles = os.listdir(path)
    for file in pklfiles:
        filepath = path + file
        phdiagram = read_pickle_object(filepath)['diagram']
        for i in range(0, len(phdiagram)):
            phdiagram[i] = np.delete(phdiagram[i], -1, axis=1)
        bary = lagrangian_barycenter(pdiagset=phdiagram) # this function may take time if phdiagram is large
        if phdiagram == []:
            print(filepath)
        else:
            diagramlist.append(bary)
    index = len(diagramlist)    
    return diagramlist, index 

diagramlist = []
namelist = ['breastmnist class',
            'chestmnist class',
            'dermamnist class',
            'octmnist class',
            'organmnist axial class',
            'organmnist coronal class',
            'organmnist sagittal class',
            'pathmnist class',
            'pneumoniamnist class',
            'retinamnist class']
indxlist = []

# obtain list od diagrams
for name in namelist:
    diagramlist, index = make_diagramlist(name, diagramlist)
    indxlist.append(index)
    
# log inf error
ep = 1e-4

# compute distance matrix
Dmatrix = np.zeros((index,index))

for i in range(0,index):
    for j in range(0,index):
        dataA = np.log(diagramlist[i]+ep)
        dataB = np.log(diagramlist[j]+ep)
        dilaInv, bottleneckDistance, minBottleneckDistance = dist.myDistance(dataA, dataB, 100)
        Dmatrix[i,j] = minBottleneckDistance

        
# mds visualization
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=10,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit(Dmatrix).embedding_

indxlist = [0] + indxlist

plt.figure()
for i in range(len(namelist)):
    plt.scatter(pos[indxlist[i]:indxlist[i+1]+1, 0], pos[indxlist[i]:indxlist[i+1]+1, 1], label = namelist[i])
plt.legend()
plt.savefig('mds_visual.png')


# umap visualization
RANDOM_SEED = 1 
reducer = umap.UMAP(
    n_neighbors = 5,       # default: 15
    n_components = 2,       # 2D atlas
    metric = 'precomputed', # we have already computed the pairwise distances
    min_dist = .05,         # default: 0.1
    spread = 1,             # default: 1
    random_state = RANDOM_SEED)

embedding = reducer.fit_transform(Dmatrix) 

# Plot the UMAP atlas 
plt.figure()
for i in range(len(namelist)):
    plt.scatter(embedding[indxlist[i]:indxlist[i+1]+1, 0], embedding[indxlist[i]:indxlist[i+1]+1, 1], label = namelist[i])
plt.legend()
plt.savefig('umap_visual.png')

