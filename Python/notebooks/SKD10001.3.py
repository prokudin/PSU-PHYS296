# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 20:21:39 2018

@author: Dolam
"""
# I may be flying south when Im supposed to be flying north.....Should I turn around???
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
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
Q2Bin=np.array([1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 5.0, 15.0])
zBin= np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.1])

# creates index values for final plot matrix
i = np.arange(81).reshape(9,9)
value = dat['value']
qT2 = dat['qT2']
delta = dat['delta']
xClas=range(len(xBin)-1)
Q2Clas=range(len(Q2Bin)-1)
zClas=range(len(zBin)-1)
ind = np.arange(336)

dat['xClas'] = pd.cut(dat['x'], xBin, labels=xClas)
dat['xBin'] = pd.cut(dat['x'], xBin,labels = False, retbins=0)
xBind = dat['xClas']
xBinned = pd.DataFrame({'xBin': dat['xBin']},index = ind)

dat['Q2Clas'] = pd.cut(dat['Q2'], Q2Bin, labels=Q2Clas)
dat['Q2Bin'] = pd.cut(dat['Q2'], Q2Bin,labels = False, retbins=0)
Q2Bind = dat['Q2Clas']
Q2Binned = pd.DataFrame({'Q2Bin' : dat['Q2Bin']}, index = ind)

dat['zClas'] = pd.cut(dat['z'], zBin, labels=zClas)
dat['zBin'] = pd.cut(dat['z'], zBin,labels = False, retbins=0)
zBind = dat['zBin']
zBinned = pd.DataFrame({'zBin' : dat['zBin']}, index = ind)
data_binvalue = dat['value'].as_matrix(columns=None)
data_binqT2 = dat['qT2'].as_matrix(columns=None)
data_bindelta = dat['delta'].as_matrix(columns=None)
data_binz = dat['z'].as_matrix(columns=None)
data_binx = dat['x'].as_matrix(columns=None)

# Creates the index needed to create f and g DataFrames
ind = np.arange(336)
ind2 = np.arange(9)
X = dat['x']

h = pd.DataFrame({'xBin': xBind, 'Q2Bin': Q2Bind, 'zBin' : zBind},index=ind)
# Creates a DataFrame with qT2 as x, value as y
g = pd.DataFrame({'x': data_binx, 'value': data_binvalue, 'z' : data_binz, 'delta' : data_bindelta, 'qT2' :
    data_binqT2},index=ind)
xbin0 = xBinned.query('xBin == 0')
x0ind  = xbin0.index.values
new_xbin0 = X.iloc[x0ind]
xbin1 = xBinned.query('xBin == 1')
x1ind  = xbin1.index.values
new_xbin1 = X.iloc[x1ind]
xbin2 = xBinned.query('xBin == 2')
x2ind  = xbin2.index.values
new_xbin2 = X.iloc[x2ind]
xbin3 = xBinned.query('xBin == 3')
x3ind  = xbin3.index.values
new_xbin3 = X.iloc[x3ind]
xbin4 = xBinned.query('xBin == 4')
x4ind  = xbin4.index.values
new_xbin4 = X.iloc[x4ind]
xbin5 = xBinned.query('xBin == 5')
x5ind  = xbin5.index.values
new_xbin5 = X.iloc[x5ind]
xbin6 = xBinned.query('xBin == 6')
x6ind  = xbin6.index.values
new_xbin6 = X.iloc[x6ind]
xbin7 = xBinned.query('xBin == 7')
x7ind  = xbin7.index.values
new_xbin7 = X.iloc[x7ind]
xbin8 = xBinned.query('xBin == 8')
x8ind  = xbin8.index.values
new_xbin8 = X.iloc[x8ind]

#Learning to combine multiple arrays back into one, This is = to ind up above
xbindatind=np.concatenate((x0ind,x1ind,x2ind,x3ind,x4ind,x5ind,x6ind,x7ind,x8ind), axis=0)

#Here is binned data back in a DataFrame idk why I did this
xbindat = pd.DataFrame({'0':new_xbin0,'1':new_xbin1,'2':new_xbin2,'3':new_xbin3,
                        '4':new_xbin4,'5':new_xbin5,'6':new_xbin6,'7':new_xbin7,
                        '8':new_xbin8},index=ind)
    
Q2 = dat['Q2']    
Q2bin0 = Q2Binned.query('Q2Bin == 0')
q0ind  = Q2bin0.index.values
new_Q2bin0 = Q2.iloc[q0ind]
Q2bin1 = Q2Binned.query('Q2Bin == 1')
q1ind  = Q2bin1.index.values
new_Q2bin1 = Q2.iloc[q1ind]
Q2bin2 = Q2Binned.query('Q2Bin == 2')
q2ind  = Q2bin2.index.values
new_Q2bin2 = Q2.iloc[q2ind]
Q2bin3 = Q2Binned.query('Q2Bin == 3')
q3ind  = Q2bin3.index.values
new_Q2bin3 = Q2.iloc[q3ind]
Q2bin4 = Q2Binned.query('Q2Bin == 4')
q4ind  = Q2bin4.index.values
new_Q2bin4 = Q2.iloc[q4ind]
Q2bin5 = Q2Binned.query('Q2Bin == 5')
q5ind  = Q2bin5.index.values
new_Q2bin5 = Q2.iloc[q5ind]
Q2bin6 = Q2Binned.query('Q2Bin == 6')
q6ind  = Q2bin6.index.values
new_Q2bin6 = Q2.iloc[q6ind]
Q2bin7 = Q2Binned.query('Q2Bin == 7')
q7ind  = Q2bin7.index.values
new_Q2bin7 = Q2.iloc[q7ind]
Q2bin8 = Q2Binned.query('Q2Bin == 8')
q8ind  = Q2bin8.index.values
new_Q2bin8 = Q2.iloc[q8ind]

Q2bindat = pd.DataFrame({'0':new_Q2bin0,'1':new_Q2bin1,'2':new_Q2bin2,'3':new_Q2bin3,
                        '4':new_Q2bin4,'5':new_Q2bin5,'6':new_Q2bin6,'7':new_Q2bin7,
                        '8':new_Q2bin8},index=ind)

Z = dat['z']
zbin0 = zBinned.query('zBin == 0')
z0ind  = zbin0.index.values
new_zbin0 = Z.iloc[z0ind]
zbin1 = zBinned.query('zBin == 1')
z1ind  = zbin1.index.values
new_zbin1 = Z.iloc[z1ind]
zbin2 = zBinned.query('zBin == 2')
z2ind  = zbin2.index.values
new_zbin2 = Z.iloc[z2ind]
zbin3 = zBinned.query('zBin == 3')
z3ind  = zbin3.index.values
new_zbin3 = Z.iloc[z3ind]
zbin4 = zBinned.query('zBin == 4')
z4ind  = zbin4.index.values
new_zbin4 = Z.iloc[z4ind]
zbin5 = zBinned.query('zBin == 5')
z5ind  = zbin5.index.values
new_zbin5 = Z.iloc[z5ind]

zbindat = pd.DataFrame({'0':new_zbin0,'1':new_zbin1,'2':new_zbin2,'3':new_zbin3,
                        '4':new_zbin4,'5':new_zbin5},index=ind)

zdatcheck = pd.DataFrame({'z': dat['z']},index=ind)
xdatcheck = pd.DataFrame({'x': dat['x']},index=ind)
# A look at the data that was able to be binned

#for index, row in xdatcheck.iterrows():
#    print(row == 1)



QT2_0 = qT2.iloc[x0ind]
value_0 = value.iloc[q0ind]
delta_0 = delta.iloc[z0ind]




