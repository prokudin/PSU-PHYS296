#!/usr/bin/env python
import sys,os
import numpy as np
import fCJLO,fCJNLO

class CJ:

  def __init__(self,conf):
    #root='./'):

    root=conf['path2CJ']
    if root.endswith('/')==False: root+='/'

    if conf['order']=='LO':
      fCJLO.root.root=root.ljust(255)
      fCJLO.setcj(400) # CJ15central
      self.fCJ=fCJLO
    elif conf['order']=='NLO': 
      fCJNLO.root.root=root.ljust(255)
      fCJNLO.setcj(500) # CJ15central
      self.fCJ=fCJNLO
    else:
      raise ValueError('ERR at %s.__init__'%self.__class__.__name__) 
    self.storage={}

  def get_f(self,x,Q2):
    """
    out: g,u,ub,d,db,s,sb,c,cb,b,bb
    """
    if (x,Q2) not in self.storage:
      self.storage[(x,Q2)]=np.array(self.fCJ.cjpdf_all(x,Q2**0.5))
    return self.storage[(x,Q2)]

if __name__=='__main__':

  conf={}
  conf['path2CJ']='./'

  conf['order']='LO'
  cjLO=CJ(conf)

  conf['order']='NLO'
  cjNLO=CJ(conf)

  x=0.1
  Q2=1.62
  print cjLO.get_f(x,Q2)
  print cjNLO.get_f(x,Q2)



