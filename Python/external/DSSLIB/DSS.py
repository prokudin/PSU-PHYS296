#!/usr/bin/env python
import sys,os
import numpy as np
import dssfpi,dssfk,dssfh

class DSS:
  
  def __init__(self,conf):

    root=conf['path2DSS']
    if root.endswith('/')==False: root+='/'

    self.io=0
    dssfpi.root.root=root.ljust(255)
    dssfk.root.root=root.ljust(255)
    dssfh.root.root=root.ljust(255)
    
    self.pi={0:{},1:{},-1:{}}
    self.k={0:{},1:{},-1:{}}
    self.h={0:{},1:{},-1:{}}

  def _get_f(self,z,Q2,storage,func,ih,ic):
    if (z,Q2) not in storage[ic]:
      storage[ic][(z,Q2)]=np.array(func(ih,ic,self.io,z,Q2))/z
    return storage[ic][(z,Q2)]

  def get_f(self,z,Q2,hadron):
    """
    out: g,u,ub,d,db,s,sb,c,cb,b,bb
    """
    if   hadron=='pi+':return self._get_f(z,Q2,self.pi,dssfpi.fdss,1,1)
    elif hadron=='pi-':return self._get_f(z,Q2,self.pi,dssfpi.fdss,1,-1)
    elif hadron=='pi0':return self._get_f(z,Q2,self.pi,dssfpi.fdss,1,0)
    elif hadron=='k+':return self._get_f(z,Q2,self.k,dssfk.fdss,2,1)
    elif hadron=='k-':return self._get_f(z,Q2,self.k,dssfk.fdss,2,-1)
    elif hadron=='k0':return self._get_f(z,Q2,self.k,dssfk.fdss,2,0)
    elif hadron=='h+':return self._get_f(z,Q2,self.h,dssfh.fdss,4,1)
    elif hadron=='h-':return self._get_f(z,Q2,self.h,dssfh.fdss,4,-1)
    elif hadron=='h0':return self._get_f(z,Q2,self.h,dssfh.fdss,4,0)

if __name__=='__main__':

  conf={}
  conf['path2DSS']='./'

  dss=DSS(conf)
  print dss.get_f(0.5,10,'pi+')
  print dss.get_f(0.5,10,'pi-')


