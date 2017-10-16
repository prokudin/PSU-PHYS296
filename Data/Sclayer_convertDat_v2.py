# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 08:33:48 2017
PSU-Phys_296: Idependent Study
Purpous: 
    Read and convert data set
    Plot data
@author: Aardvark
"""
## Reading data as .dat 
from numpy import genfromtxt as gft # Import package
from pandas import DataFrame as df # Import package

L=open('clas_data.dat').readlines()     # read each line
L=[l.strip() for l in L]                # remove white space between values in each line
H=L[1].replace('#PhT','PhT').split()    # save the second line as H and delete the hashtag 
dat = gft('clas_data.dat')              # Read file from working dir as numpi array
data = df(dat, columns = H)             # save data as a panda data frame with headers

# saving as csv
data.to_csv("class_data_2.csv",index=False)    # save as csv with out index vaues


