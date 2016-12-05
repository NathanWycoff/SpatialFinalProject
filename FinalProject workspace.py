# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 18:40:18 2016

@author: nathan
"""

import numpy as np
import matplotlib.pyplot as plt

####
###Data Generation steps
####
#Define our difEQ's (with an analytic solution)
dydt = lambda y,t: y
y = lambda t: np.exp(t)

#Generate some noise corrupted data
n = 10
t_o = np.linspace(0,2,num=n)
y_o = np.array([y(x) for x in t_o]) + np.random.normal(size=n)

#Plot the data
plt.plot(t_o, y_o)

####
#Empirical Covariance function
####
#Params
bin_count = 5#How many distance bins for the covariogram?

#Create distance matrix
dist_mat = np.array([[np.linalg.norm(x1-x2) for x1 in t_o] for x2 in t_o])

#Get bins
d_extrema = np.percentile(dist_mat.flatten(), [1,100])
d_range = d_extrema[1] - d_extrema[0]
bins = [d_extrema[0] + d_range/bin_count * (x+1) for x in range(bin_count)]




####
###Ordinary GP with gaussian kernel function
####
k = lambda x1, x2: np.exp(-np.linalg.norm(x1 - x2))

#Build covariance matrix
K = np.array([[k(x1,x2) for x1 in t_o] for x2 in t_o])

#Predict for t = 0.5
k = np.array([k(0.5,x) for x in t_o])

#Prediction 
L = np.linalg.cholesky(K)
alpha = np.linalg.solve(L.T, np.linalg.solve(L,y_o))
pred = np.dot(k.T, alpha)