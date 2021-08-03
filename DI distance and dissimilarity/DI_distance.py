# -*- coding: utf-8 -*-
"""
Created on Tue May 25 07:13:44 2021

@author: Yueqi
"""

import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt

# main function

def myDistance(X, Y, num_partitions):
    
    ## character point
    indX = np.argmax(X[:,1]-X[:,0])
    indY = np.argmax(Y[:,1]-Y[:,0])
    chaX = X[indX,:]
    chaY = Y[indY,:]
    
    ## define the searching interval
    B = 1e3 ## translation parameter. This is a bug in gudhi that they donot accept negative inputs
    oriDist = gd.bottleneck_distance(X+B,Y+B)
    print("original bottleneck distance", oriDist)
    upBound = max(Y[:,1])-chaX[1]+oriDist
    lowBound = chaY[1]-max(X[:,1])-oriDist
    
    ## define shifts
    shifts = np.linspace(lowBound, upBound, num_partitions)
    shiftInv = np.append(shifts, 0)
    bottleneckDistance = []
    
    ## compute bottleneck distance
    for scale in shiftInv:
        shiftX = scale+X
        dist = gd.bottleneck_distance(shiftX+B, Y+B)
        bottleneckDistance.append(dist)

    minBottleneckDistance = min(bottleneckDistance)
    print("dilation-invariant bottleneck distance (using log) is", minBottleneckDistance)
    
    return shiftInv, bottleneckDistance, minBottleneckDistance


# load data

##file_list = ['activity_net_euclidean_cosine.npz','activity_net_euclidean.npz']
##file_list = ['activity_net_euclidean_expmap_poincare.npz','activity_net_poincare.npz']
#file_list = ['mammals_euclidean_cosine.npz','mammals_euclidean.npz']
##file_list = ['mammals_euclidean_expmap_poincare.npz','mammals_poincare.npz']
#
#dist_list = [np.load(u)['arr_0'] for u in file_list]
#
#X = dist_list[0][0]
#Y = dist_list[1][0]
#
#X = np.delete(X, -1, axis = 1)
#Y = np.delete(Y, -1, axis = 1)
#
### avoid -inf in log map
#ep = 1e-10
#X = np.log(X+ep)
#Y = np.log(Y+ep)
#
#shiftInv, bottleneckDistance, minBottleneckDistance = myDistance(X, Y, 100)
#
#plt.figure()
#label = ['different shifts','no shift','optimal shift']
#plt.xlabel("log shift parameters")
#plt.ylabel("bottleneck distances")
#plt.plot(shiftInv[0:len(shiftInv)-1], bottleneckDistance[0:len(bottleneckDistance)-1])
## the last element corresponds to no dilation
#plt.axhline(y=bottleneckDistance[-1], color='r', linestyle='-.')
#plt.axhline(y=minBottleneckDistance, color='g', linestyle='-')
#plt.legend(label, loc='best')