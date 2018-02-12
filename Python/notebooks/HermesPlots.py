# -*- coding: utf-8 -*-
"""
Spyder Editor

This file was made for ploting Hermes data.
"""
## Packages
import os
import pandas as pd
import numpy as np
from pandas import read_excel as read
from pandas import DataFrame as df
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

## Set working directory (Windows)
os.chdir('../') # move dir back one 
current = os.getcwd() # get current dir
print "Moved Working directory back  %s" % current # check current dir

path = current+'\\database\\sidis\\expdata\\' #for Windows

#path = current+'/database/sidis/expdata/' # for Linux

os.chdir(path)
# Check 
newDir = os.getcwd() # get current dir
print "Working directory for data %s" % newDir # check current dir

## Grabing list of files
data=os.listdir('./') # list all files in dir
data=[files for files in data if files.endswith('.xlsx') and files.startswith('1')] # list of COMPASS data
print "Data files retrieved %s" % data


## Reading and resructuring data

dat = df(read(data[0]))

# Calculate 
dat["delta"] = np.sqrt(dat["stat_u"]**2.0) # measurment error
dat["qT"] = dat["pT"]/dat["z"]
dat["qT2"] = dat["qT"]**2 #

##Binning data
xBin=[0.023,0.04,0.055,0.075,0.1,0.14,0.2,0.3,0.4,0.6]
Q2Bin=[1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 5.0, 15.0]
zBin= [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.1]

xClas=range(len(xBin)-1)
Q2Clas=range(len(Q2Bin)-1)
zClas=range(len(zBin)-1)

dat['xClas'] = pd.cut(dat['x'], xBin, labels=xClas)
dat['xBin'] = pd.cut(dat['x'], xBin)

dat['Q2Clas'] = pd.cut(dat['Q2'], Q2Bin, labels=Q2Clas)
dat['Q2Bin'] = pd.cut(dat['Q2'], Q2Bin)

dat['zClas'] = pd.cut(dat['z'], zBin, labels=zClas)
dat['zBin'] = pd.cut(dat['z'], zBin)

## Seting subplot perameters
Zcolor=["red","green","blue","orange",'black','purple','pink']
Zmark=['o', 'v', '^', '<', '>','*','X']
Zline=['--','--','--','--','--','--','--']

shairYax=[[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[7,1],[8,1],[9,1]] #Subplots with y-axes ticks
shairXax=[[9,1],[9,2],[9,3],[9,4],[9,5],[9,6],[9,7],[9,8],[9,9]] #Subplots with x-axes ticks

# Trouble of many kind, some come from above, some come from behine,

k=0
hermes={}
hermesZ={}
for n in xClas:
    for m in Q2Clas:
        data_bin=dat.('xClas=="%s" and Q2Clas=="%s"' %(n,m),inplace=False)
        if data_bin != []:
            hermes={"%s%s%s" %(k,n,m): data_bin}
        print data_bin.head()
        print (k,n,m)
        for i in zClas:
            data_Z = [data_bin.zClas==i]
            if data_Z != []:
                hermesZ={"%s%s%s%s" %(k,n,m,i): data_Z}
        k+=1
        
## But I bought a big batt, you'll see, my trouble are gonna have troubles with me. 

## The plot
fig=plt.figure(figsize=(15, 15),facecolor="gray") # figsize=wxh in inches
subMat = gridspec.GridSpec(ncols=9, nrows=9, wspace=0.0, hspace=0.0)

ax=fig.add_axes([0,0,1,1])
#ax.yaxis.set_ticks([ 0.,0.1245,0.2755,0.427,0.57752,0.729,0.8805])
#ax.xaxis.set_ticks([0,0.1285,0.2255,0.322,0.419,0.516,0.613,0.709,0.71+0.097,0.9])
#ax.set_yticklabels([0]+Q2Bin)
#ax.set_xticklabels([0]+xBin)

# Set title and axis labels
#ax.title("COMPASS Data")
#ax.set_ylabel(r"$Q^2$ bins",rotation="horizontal")
#ax.set_xlabel(r"$x$ bins")

# Set legend
#Zpatch1 = mpatches.Patch(color=Zcolor[0], label='z=0.2')
#Zpatch2 = mpatches.Patch(color=Zcolor[1], label='z=0.3')
#Zpatch3 = mpatches.Patch(color=Zcolor[2], label='z=0.4')
#Zpatch4 = mpatches.Patch(color=Zcolor[3], label='z=0.6')
#ax.legend(handles=[Zpatch1,Zpatch2,Zpatch3,Zpatch4],loc='upper left')
#ax.grid()

for n,N in zip(range(len(xBin)),xClas): 
    for m,M in zip(range(len(Q2Bin)),Q2Clas):
        data_bin=dat.query('xClas=="%s" and Q2Clas=="%s"' %(N,M))
        k=0 #counter
        for i in zBin:
            ax = fig.add_subplot(subMat[int(n),int(m)])
            ax.set_yscale("log")
            ax.errorbar(data_bin.qT2[data_bin.zBin==i],data_bin.value[data_bin.zBin==i],#x,y
                        data_bin.delta[data_bin.zBin==i],#errorbars
                        color=Zcolor[k],marker=Zmark[k],linestyle=Zline[k],linewidth=0,markersize=2)#line properties
            if [int(n),int(m)] not in shairYax:
                ax.set_yticklabels('')
            if [int(n),int(m)] not in shairXax:
                ax.set_xticklabels('')
            k+=1 #add one to counter
            ax.set_xlabel(r"$q_T^2$")
            
            
            
            
            
            
            
            
            
            
            
            
            
            