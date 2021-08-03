# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 09:37:45 2021

@author: Administrator
"""

import numpy as np
import DI_dissimilarity as dism

# load data
    
file_list = ['activity_net_euclidean_cosine.npz','activity_net_euclidean.npz',
             'activity_net_euclidean_expmap_poincare.npz','activity_net_poincare.npz',
             'mammals_euclidean_cosine.npz','mammals_euclidean.npz',
             'mammals_euclidean_expmap_poincare.npz','mammals_poincare.npz']

dist_list = [np.load(u)['arr_0'] for u in file_list]

# asymmetric version

ADsim = np.zeros((8,8))
for i in range(0,8):
    for j in range(0,8):
        dataA = np.delete(dist_list[i][0], -1, axis = 1)
        dataB = np.delete(dist_list[j][0], -1, axis = 1)
        dilaInv, bottleneckDistance, minBottleneckDistance = dism.myDissimilarity(dataA, dataB, 10)
        ADsim[i,j] = minBottleneckDistance
        
# symmetric version

Dsim = np.zeros((8,8))
for i in range(0,8):
    for j in range(0,8):
        dataA = np.delete(dist_list[i][0], -1, axis = 1)
        dataB = np.delete(dist_list[j][0], -1, axis = 1)
        dilaInv, bottleneckDistance, minBottleneckDistance = dism.myDissimilarity(dataA, dataB, 10)
        rdilaInv, rbottleneckDistance, rminBottleneckDistance = dism.myDissimilarity(dataB, dataA, 10)
        Dsim[i,j] = 0.5*(minBottleneckDistance + rminBottleneckDistance)

