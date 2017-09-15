# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 08:33:48 2017
PSU-Phys_296: Idependent Study
Purpous: 
    Read and convert data set
    Plot data
@author: Aardvark
"""
# Reading data as .dat 
from numpy import genfromtxt as gft # Import package

dat = gft('clas_data.dat')          # Read file from working dir

# saving as csv
from pandas import DataFrame as df # Import package

data = df(dat)                              # data as a panda data frame
data.to_csv("clas_data.csv",index=False)    # save as csv with out index vaues

# plot Pht with random and systimatic error
import matplotlib.pyplot as plt # Import package

plt.figure() # this seems to be usless
plt.errorbar(dat[:,0],dat[:,1],dat[:,2],dat[:,3],fmt='o') # bild scatter with error bars
plt.title('Figure 1') # add title
plt.show() # seems to be usless
plt.clf() # should be clearing the figure but it is only ploting in my console so, usless?

# plot Pht with random error

plt.errorbar(dat[:,0],dat[:,1],yerr=dat[:,2],fmt='o',capsize=15)
plt.title('Figure 2')
