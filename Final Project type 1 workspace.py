# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 17:13:58 2016

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

#Generate some noise corrupted data, this is what we want to do analysis on.
n = 10
t_o = np.linspace(0,2,num=n)
y_o = np.array([y(x) for x in t_o]) + np.random.normal(size=n)

#Plot the data
plt.plot(t_o, y_o)