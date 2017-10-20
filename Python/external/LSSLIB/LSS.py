#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import fLSS

class LSS:

  def __init__(self,conf):

    root=conf['path2LSS']
    if root.endswith('/')==False: root+='/'
    fLSS.root.root=root.ljust(255)
    self.iset=1 # positive glue
    self.storage={}

  def _get_f(self,x,Q2):
    """
    out: g,u,ub,d,db,s,sb,c,cb,b,bb
    """
    a=np.array(fLSS.lss2014(self.iset,x,Q2))
    g=a[3]
    u=a[0]-a[2]
    d=a[1]-a[2]
    s=a[2]/2
    c=0.0  #for now
    b=0.0  #for now
    ub=s   #for now
    db=s   #for now
    sb=s   #for now
    cb=0.0 #for now
    bb=0.0 #for now
    return np.array([g,u,ub,d,db,s,sb,c,cb,b,bb])
 
  def get_f(self,x,Q2):
    if (x,Q2) not in self.storage:
      self.storage[(x,Q2)]=self._get_f(x,Q2)
    return self.storage[(x,Q2)]

if __name__=='__main__':

  conf={}
  conf['path2LSS']='./'
  lss=LSS(conf)
  x=0.1
  Q2=1.62
  print lss.get_f(x,Q2)  




