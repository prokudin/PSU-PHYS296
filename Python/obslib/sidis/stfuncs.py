#!/usr/bin/env python
import sys,os
import numpy as np
import math
from tools.tools import load_config
from external.CJLIB.CJ import CJ
from external.LSSLIB.LSS import LSS
from external.DSSLIB.DSS import DSS
from qcdlib.tmdlib import PDF,PPDF,FF
from qcdlib.aux import AUX
import matplotlib.pyplot as plt

class STFUNCS:

  def __init__(self,conf):
    self.aux=conf['aux']  
    self.conf=conf
    eu2,ed2=4/9.,1/9. 
    self.e2=[]
    self.e2.append(0)   # g
    self.e2.append(eu2) # u
    self.e2.append(eu2) # ub
    self.e2.append(ed2) # d
    self.e2.append(ed2) # db
    self.e2.append(ed2) # s
    self.e2.append(ed2) # sb
    self.e2.append(0)   # c
    self.e2.append(0)   # cb
    self.e2.append(0)   # b
    self.e2.append(0)   # bb
    self.e2=np.array(self.e2)

    self.Mh={}
    self.Mh['pi+']=self.aux.Mpi
    self.Mh['pi-']=self.aux.Mpi
    self.Mh['k+']=self.aux.Mk
    self.Mh['k-']=self.aux.Mk

    self.D={}
    self.D[1] ={'k1':'pdf','k2':'ff'}
    self.D[2] ={'k1':'ppdf','k2':'ff'}


  def get_K(self,i,x,Q2,z,pT,wq,k1,k2,target,hadron):
    if   i==1: return x
    elif i==2: return x
    


  def get_wq(self,z,k1,k2,target,hadron):
    return z**2*np.abs(self.conf[k1].widths[target]) + np.abs(self.conf[k2].widths[hadron])

  def get_gauss(self,z,pT,wq):
    return np.exp(-pT**2/wq)/(np.pi*wq)

  def get_FX(self,i,x,z,Q2,pT,target,hadron):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    mu2=Q2
    F=self.conf[k1].get_C(x,mu2,target)
    D=self.conf[k2].get_C(z,mu2,hadron)
    wq=self.get_wq(z,k1,k2,target,hadron)
    gauss=self.get_gauss(z,pT,wq) 
    K=self.get_K(i,x,Q2,z,pT,wq,k1,k2,target,hadron)
    return np.sum(self.e2*K*F*D*gauss)

  def FLL(self,x,Q2,y,z,pT,target,hadron):
    coupling=1.0/137
    p2=y*(1.0-y/2)/(1.0-y+y**2/2)
    factor=coupling**2/(x*y*Q2)*(1-y+y**2/2)
    return factor*(p2*self.get_FX(2,x,z,Q2,pT,target,hadron))

  def FUU(self,x,Q2,y,z,pT,target,hadron):
    coupling=1.0/137
    factor=coupling**2/(x*y*Q2)*(1-y+y**2)
    return factor*(self.get_FX(1,x,z,Q2,pT,target,hadron))

if __name__=='__main__':

  conf={}
  conf['path2CJ'] ='../../external/CJLIB'
  conf['path2LSS']='../../external/LSSLIB'
  conf['path2DSS']='../../external/DSSLIB'

  conf['order']='LO'
  conf['aux']  =AUX()
  conf['_pdf'] =CJ(conf)
  conf['_ppdf']=LSS(conf)
  conf['_ff']  =DSS(conf)

  conf['pdf']=PDF(conf)
  conf['ppdf']=PPDF(conf)
  conf['ff']=FF(conf)


  stfuncs=STFUNCS(conf)
  x=0.25
  z=0.5
  Q2=2.4
  mu2=2.0
  E=11.0
  m=0.938
  pT=0.3
  y=Q2/(2*m*E*x)
  target='p'
  hadron='pi+' 
  for i in range(1,2): print i,stfuncs.get_FX(i,x,z,Q2,pT,target,hadron)
  print stfuncs.FLL(x,Q2,y,z,pT,target,hadron)
  print stfuncs.FUU(x,Q2,y,z,pT,target,hadron)
  for j in range(1,36): plt.plot([j/37.] , [stfuncs.FLL(j/37.,Q2,y,z,pT,target,hadron)], 'ro')

  plt.xlabel('x')
  plt.ylabel('FLL')
  plt.title('FLL strucrure function at  $Q^2=3.6$')

  plt.show()
  
