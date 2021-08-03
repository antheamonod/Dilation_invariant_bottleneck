# -*- coding: utf-8 -*-
"""
Created on Sun May  9 14:48:21 2021

@author: Yueqi
"""

import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt


# main function
def myDissimilarity(X, Y, num_partitions):
    
    perX = max(X[:,1]-X[:,0])
    indX = np.argmax(X[:,1]-X[:,0])
    perY = max(Y[:,1]-Y[:,0])
    indY = np.argmax(Y[:,1]-Y[:,0])
        
    ## character point
    chaX = X[indX,:]
    chaY = Y[indY,:]
    
    ## define the searching interval
    oriDist = gd.bottleneck_distance(X,Y)
    print("original bottleneck distance", oriDist)
    upBound = min((perY+2*oriDist)/perX, (max(Y[:,1])+oriDist)/chaX[1])
    lowBound = max((perY-2*oriDist)/perX, (chaY[1]-oriDist)/max(X[:,1]))
    
    ## define partitions
    dilations = np.linspace(lowBound, upBound, num_partitions)
    dilaInv = np.append(dilations, 1)
    
    bottleneckDistance = []
    
    ## compute bottleneck distance
    for scale in dilaInv:
        dilaX = scale*X
        dist = gd.bottleneck_distance(dilaX, Y)
        bottleneckDistance.append(dist)
        
    minBottleneckDistance = min(bottleneckDistance)
    #loc = np.argmin(bottleneckDistance)
    #optdila = dilaInv[loc]
    print("dilation-invariant dissimilarity (left-to-right) is", minBottleneckDistance)
    
    return dilaInv, bottleneckDistance, minBottleneckDistance


# load data
    
#file_list = ['activity_net_euclidean_cosine.npz','activity_net_euclidean.npz']
#file_list = ['activity_net_euclidean_expmap_poincare.npz','activity_net_poincare.npz']
#file_list = ['mammals_euclidean_cosine.npz','mammals_euclidean.npz']
#file_list = ['mammals_euclidean_expmap_poincare.npz','mammals_poincare.npz']
#
#dist_list = [np.load(u)['arr_0'] for u in file_list]
#
#X = dist_list[0][0]
#Y = dist_list[1][0]
#
#X = np.delete(X, -1, axis = 1)
#Y = np.delete(Y, -1, axis = 1)
#
#dilaInv, bottleneckDistance, minBottleneckDistance = myDissimilarity(X, Y, 100)
#
#plt.figure(figsize=(10,8))
#label = ['Different Dilations','No Dilation','Optimal Dilation']
#plt.rcParams.update({"font.size": 20, "font.family":"times new roman"})
#plt.xlabel("Dilation Parameters")
#plt.ylabel("Bottleneck Distances")
#plt.plot(dilaInv[0:len(dilaInv)-1], bottleneckDistance[0:len(bottleneckDistance)-1])
## the last element corresponds to no dilation
#plt.axhline(y=bottleneckDistance[-1], color='r', linestyle='-.')
#plt.axhline(y=minBottleneckDistance, color='g', linestyle='-')
#plt.legend(label, loc='best')
