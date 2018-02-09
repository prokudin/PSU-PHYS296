#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import fGRSV

class GRSV:

  def __init__(self,conf):

    root=conf['path2GRSV']
    if root.endswith('/')==False: root+='/'
    fGRSV.root.root=root.ljust(255)
    self.iset=3 # LO standard
    self.storage={}

  def _get_f(self,x,Q2):
    """
    out: U, D, UB, DB, ST, GL, G1P, G1N
    """
    a=np.array(fGRSV.parpol(self.iset,x,Q2))
    g=a[4]
    u=a[0]
    d=a[1]
    s=a[3]
    c=0.0   
    b=0.0   
    ub=a[2]    
    db=a[3]   
    sb=s    
    cb=0.0  
    bb=0.0 
    return np.array([g,u,ub,d,db,s,sb,c,cb,b,bb])/x
 
  def get_f(self,x,Q2):
    if (x,Q2) not in self.storage:
      self.storage[(x,Q2)]=self._get_f(x,Q2)
    return self.storage[(x,Q2)]

if __name__=='__main__':

  conf={}
  conf['path2GRSV']='./'
  grsv=GRSV(conf)
  x=0.1
  Q2=1.62
  print grsv.get_f(x,Q2)  




