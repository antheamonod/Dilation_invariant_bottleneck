# -*- coding: utf-8 -*-
"""
Created on Tue May 25 07:21:19 2021

@author: Yueqi
"""

import numpy as np
import main_algorithm as mydist
import gudhi as gd

#file_list = ['activity_net_euclidean_cosine.npz','activity_net_euclidean.npz']
#file_list = ['activity_net_euclidean_expmap_poincare.npz','activity_net_poincare.npz']
file_list = ['mammals_euclidean_cosine.npz','mammals_euclidean.npz']
#file_list = ['mammals_euclidean_expmap_poincare.npz','mammals_poincare.npz']

# Read files
dist_list =[np.load(u)['arr_0'] for u in file_list]

X = dist_list[0][0]
Y = dist_list[1][0]

X = np.delete(X, -1, axis = 1)
Y = np.delete(Y, -1, axis = 1)

# subsamples. the whole data will take 1 week
X = X[:32,:]
Y = Y[:32,:]

# log map
ep=1e-10
X = np.log(X+ep)
Y = np.log(Y+ep)

# translation
B = 1e3
ori = gd.bottleneck_distance(X+B,Y+B)
print('original distance is ', ori)

X = [tuple(x) for x in X.tolist()]
Y = [tuple(x) for x in Y.tolist()]

# main function
DShift = mydist.shifted_bottleneck_distance(X,Y, analysis = True)

print('Don distance is ', DShift)