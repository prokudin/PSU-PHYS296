#!/usr/bin/env python
import sys,os;
#sys.path.insert(1,'../') 
from aux import AUX
import numpy as np

class ALPHAS:

  def __init__(self,conf):
    self.aux=conf['aux']

    if   conf['order']=='LO':  self.order=0
    elif conf['order']=='NLO': self.order=1
    self.conf=conf

    self.Q20 = conf['Q20']
    self.aZ  = self.aux.alphaSMZ/(4*np.pi)
    self.setup()

  def setup(self):

    self.beta=np.zeros((7,3))
    for Nf in range(3,7): 
      self.beta[Nf,0]=11.0-2.0/3.0*Nf 
      self.beta[Nf,1]=102.-38.0/3.0*Nf 
      self.beta[Nf,2]=2857.0/2.0-5033.0/18.0*Nf+325.0/54.0*Nf**2 
    
    if self.conf['alphaSmode']=='backward':
      # uses alphaS(mZ)--> backwards evolution
      self.ab=self.evolve_a(self.aux.mZ2,self.aZ,self.aux.mb2,5)
      self.ac=self.evolve_a(self.aux.mb2,self.ab,self.aux.mc2,4)
      self.a0=self.evolve_a(self.aux.mc2,self.ac,self.Q20,3)

    elif self.conf['alphaSmode']=='forward':
      self.a0=self.conf['alphaS0']/(4*np.pi)
      self.ac=self.evolve_a(self.Q20,self.a0,self.aux.mc2,3)
      self.ab=self.evolve_a(self.aux.mc2,self.ac,self.aux.mb2,4)

    # we will store all Q2 values of alphaS 
    self.storage={}

  def get_Nf(self,Q2):
    Nf=3
    if Q2>=(4.*self.aux.mc2): Nf+=1
    if Q2>=(4.*self.aux.mb2): Nf+=1
    return Nf

  def beta_func(self,a,Nf):
    betaf = -self.beta[Nf,0]
    if self.order>=1: betaf+=-a*self.beta[Nf,1]
    if self.order>=2: betaf+=-a*self.beta[Nf,2]
    return betaf*a**2

  def evolve_a(self,Q20,a,Q2,Nf):
    # Runge-Kutta implemented in pegasus  
    LR = np.log(Q2/Q20)/20.0
    for k in range(20):
      XK0 = LR * self.beta_func(a,Nf)
      XK1 = LR * self.beta_func(a + 0.5 * XK0,Nf)
      XK2 = LR * self.beta_func(a + 0.5 * XK1,Nf)
      XK3 = LR * self.beta_func(a + XK2,Nf)
      a+= (XK0 + 2.* XK1 + 2.* XK2 + XK3) * 0.166666666666666
    return a

  def get_a(self,Q2):

    if Q2 not in self.storage:
      if self.aux.mb2<=Q2:
        Q20,a0,Nf=self.aux.mb2,self.ab,5
      elif self.aux.mc2<=Q2 and Q2<self.aux.mb2: 
        Q20,a0,Nf=self.aux.mc2,self.ac,4
      elif Q2<self.aux.mc2:
        Q20,a0,Nf=self.Q20,self.a0,3
      self.storage[Q2]=self.evolve_a(Q20,a0,Q2,Nf)
    return self.storage[Q2]

  def get_alphaS(self,Q2):
    return self.get_a(Q2)*4*np.pi

if __name__=='__main__':

  conf={}
  conf['alphaSmode']='backward'
  conf['alphaS0']=0.3
  conf['mode']='truncated'
  conf['order']='NLO'
  conf['scheme']='ZMVFS'
  conf['Q20'] = 1.0
  conf['aux']=AUX()
  aS=ALPHAS(conf)

  print '========================'
  print 'test alphaS evolution'
  print '========================'
  print 'Q2=1           alphaS=%0.5f'%aS.get_alphaS(1.0)
  print 'Q2=(1+mc2)/2   alphaS=%0.5f'%aS.get_alphaS(0.5*(1.0+aS.aux.mc2))
  print 'Q2=mc2         alphaS=%0.5f'%aS.get_alphaS(aS.aux.mc2)
  print 'Q2=(mc2+mb2)/2 alphaS=%0.5f'%aS.get_alphaS(0.5*(aS.aux.mc2+aS.aux.mb2))
  print 'Q2=mb2         alphaS=%0.5f'%aS.get_alphaS(aS.aux.mb2)
  print 'Q2=(mb2+mZ2)/2 alphaS=%0.5f'%aS.get_alphaS(0.5*(aS.aux.mb2+aS.aux.mZ2))
  print 'Q2=mZ2         alphaS=%0.5f'%aS.get_alphaS(aS.aux.mZ2)
  
































