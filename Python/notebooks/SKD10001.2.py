# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 19:55:24 2018

@author: Dolam
"""

import pandas as pd
import numpy as np
from pandas import read_excel as read
from pandas import DataFrame as df
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import string
# enable plots in the notebook
#%matplotlib inline
#cc = pd.read_excel('CustomColors.xlsx')
dat=pd.read_excel('C:/Users/Dolam/Documents/Scott/1000.xlsx');

# Calculate 
dat["delta"] = np.sqrt(dat["stat_u"]**2.0) #measurment error
dat["qT"] = dat["pT"]/dat["z"]
dat["qT2"] = dat["qT"]**2

##Binning data Tick marks for overall 9x9 matrix
xBin=np.array([0.023,0.04,0.055,0.075,0.1,0.14,0.2,0.3,0.4,0.6]) 
yBin=np.array([1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 5.0, 15.0])
zBin= np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.1])

# creates index values for final plot matrix
i = np.arange(81).reshape(9,9)
ydat = dat['value']
xdat = dat['qT2']
zdat = dat['delta']
#yClas=(range(0,len(yBin))-1)
#xClas=(range(0,len(xBin))-1)
#zClas=(range(0,len(zBin))-1)

data_biny = dat['value'].as_matrix(columns=None)
data_binx = dat['qT2'].as_matrix(columns=None)
data_binz = dat['delta'].as_matrix(columns=None)

# Creates the index needed to create f and g DataFrames
ind = np.arange(336)
ind2 = np.arange(9)

# Creates a DataFrame with qT2 as x, value as y, and delta as z

h = pd.DataFrame({'x': data_binx, 'y': data_biny},index=ind)
# Creates a DataFrame with qT2 as x, value as y
g = pd.DataFrame({'x': data_binx, 'y': data_biny},index=ind)

# Same cut applied to  g
g['y'] = pd.cut(ydat, yBin, labels=None,retbins=0)
g['x'] = pd.cut(xdat, xBin, labels=None,retbins=0)
g2 = g.dropna()
print(g2)
gind = g2.index.values
newg = h.iloc[gind]
Data4plot = newg.sort_values(by=['x'])
print(Data4plot)