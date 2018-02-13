#!/usr/bin/env python
import sys,os
import numpy as np
import math
sys.path.append('./../../') #this appends/searches for the contents of the 
                  #the dir/folders that are two directories back
from tools.tools import load_config
from external.PDF.CT10 import CT10
#from external.CJLIB.CJ import CJ
from external.LSSLIB.LSS import LSS
from external.DSSLIB.DSS import DSS
from qcdlib.tmdlib import PDF,PPDF,FF,GK
from qcdlib.aux import AUX
import pylab as py
#import matplotlib.pyplot as plt
from tools.hankel.inverters import Ogata, Quad
 
class STFUNCS:  # creating a class of 

  def __init__(self,conf):
    self.aux=conf['aux']  
    self.conf=conf
    eu2,ed2=4/9.,1/9. 
    self.e2=[]   #open list
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
    self.e2=np.array(self.e2)  #create an array of charges

    self.Mh={}
    self.Mh['h+']=self.aux.Mpi
    self.Mh['h-']=self.aux.Mpi
    self.Mh['pi+']=self.aux.Mpi
    self.Mh['pi-']=self.aux.Mpi
    self.Mh['k+']=self.aux.Mk
    self.Mh['k-']=self.aux.Mk

    self.D={}   #creat a dictionalary
    self.D[1] ={'k1':'pdf','k2':'ff'}
    self.D[2] ={'k1':'ppdf','k2':'ff'}
    self.ogata = Ogata()
    self.quad = Quad()

  def get_K(self,i,x,Q2,z,pT,wq,k1,k2,target,hadron):
    if   i==1: return x
    elif i==2: return x
    

# defining functions
  def get_wq(self,z,k1,k2,target,hadron):
    return z**2*np.abs(self.conf[k1].widths[target]) + np.abs(self.conf[k2].widths[hadron])

  def get_wq_evolution(self,z,x,Q):
    wq=np.ones(len(self.e2))
    return wq * 4*z**2 * self.conf['gk'].g2 * np.log( (z*Q)/(x*self.conf['gk'].Q0) )
#    return wq * 4*z**2 * self.conf['gk'].g2 * np.log( (Q)/(self.conf['gk'].Q0) )

  def get_gauss(self,z,pT,wq):
    return np.exp(-pT**2/wq)/(np.pi*wq)

# Structure functions in pT space
  def get_FX(self,i,x,z,Q2,pT,target,hadron):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    mu2=Q2
    Q = np.sqrt(Q2)
    F=self.conf[k1].get_C(x,mu2,target)
    D=self.conf[k2].get_C(z,mu2,hadron)
    wq=self.get_wq(z,k1,k2,target,hadron)  +  self.get_wq_evolution(z,x,Q)
    gauss=self.get_gauss(z,pT,wq) 
    K=self.get_K(i,x,Q2,z,pT,wq,k1,k2,target,hadron)
    return np.sum(self.e2*K*F*D*gauss)  #sums up the contributions

  def FLL(self,x,Q2,y,z,pT,target,hadron):
    coupling=1.0/137
    p2=y*(1.0-y/2)/(1.0-y+y**2/2)
    factor=coupling**2/(x*y*Q2)*(1-y+y**2/2)
    return factor*(p2*self.get_FX(2,x,z,Q2,pT,target,hadron))

  def FUU(self,x,Q2,y,z,pT,target,hadron):
    coupling=1.0/137
    factor=coupling**2/(x*y*Q2)*(1-y+y**2)
    return factor*(self.get_FX(1,x,z,Q2,pT,target,hadron))

    
# bstar
  def bstar(self,b):
    return b/np.sqrt(1+b**2/self.conf['gk'].bmax**2)

# mub
  def mub(self,b):
    return 2*np.exp(-np.euler_gamma)/self.bstar(b)

# gk function
  def get_gk(self,b,z,x,Q):
    wq=np.ones(len(self.e2))
    return wq * self.conf['gk'].g2 * np.log(b/self.bstar(b)) * np.log( (z*Q)/(x*self.conf['gk'].Q0) )

# intrinsic widths
  def get_width(self,b,z,k1,k2,target,hadron):
    return np.abs(self.conf[k1].widths[target])/4 + np.abs(self.conf[k2].widths[hadron])/(4*z**2)
      
# Structure functions in b space
  def get_FX_b(self,i,x,z,Q2,b,target,hadron):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    mu2=self.mub(b)**2
    Q = np.sqrt(Q2)
    F=self.conf[k1].get_C(x,mu2,target)
    D=self.conf[k2].get_C(z,mu2,hadron)
    width=self.get_width(b,z,k1,k2,target,hadron)*b**2  +  self.get_gk(b,z,x,Q)
    K=self.get_K(i,x,Q2,z,pT,width,k1,k2,target,hadron)
    return 2*np.pi*np.sum(self.e2*K*F*D*np.exp(-width))  #sums up the contributions
    
# Structure function FUU in b space
  def FUU_b(self,x,Q2,y,z,b,target,hadron):
    coupling=1.0/137
    factor=coupling**2/(x*y*Q2)*(1-y+y**2)
    return factor*(self.get_FX(1,x,z,Q2,b,target,hadron))

  def FUU_q(self,x,Q2,y,z,q,target,hadron,bmin,bmax,Nmax = 10):
    nu = 0
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,b,target,hadron))
    return self.ogata.adog(w, q, nu, Nmax, bmin, bmax, lib = None)

  def FUU_q_quad(self,x,Q2,y,z,q,target,hadron,eps = 1e-3):
    nu = 0
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,b,target,hadron))
    return self.quad.quadinv(w, q, nu, eps)

#if __name__=='__main__':
#
#  conf={}
##  conf['path2CJ'] ='../../external/CJLIB'
#  conf['path2CT10'] ='../../external/PDF'
#  conf['path2LSS']='../../external/LSSLIB'
#  conf['path2DSS']='../../external/DSSLIB'
#
#  conf['order']='LO'
#  conf['aux']  =AUX()
##  conf['_pdf'] =CJ(conf)
#  conf['_pdf'] =CT10(conf)
#  conf['_ppdf']=LSS(conf)
#  conf['_ff']  =DSS(conf)
#
#  conf['pdf']=PDF(conf)
#  conf['ppdf']=PPDF(conf)
#  conf['ff']=FF(conf)
#  conf['gk']=GK(conf)
#
#
#  stfuncs=STFUNCS(conf)
#  x=0.25
#  z=0.5
#  Q2=2.4
#  mu2=2.0
#  E=11.0
#  m=0.938
#  pT=0.3
#  y=Q2/(2*m*E*x)
#  target='p'
#  hadron='pi+' 
#  for i in range(1,2): print i,stfuncs.get_FX(i,x,z,Q2,pT,target,hadron)
#  print stfuncs.FLL(x,Q2,y,z,pT,target,hadron)
#  print stfuncs.FUU(x,Q2,y,z,pT,target,hadron)
#  for j in range(1,36): plt.plot([j/37.] , [stfuncs.FLL(j/37.,Q2,y,z,pT,target,hadron)], 'ro')
#
#  plt.xlabel('x')
#  plt.ylabel('FLL')
#  plt.title('FLL strucrure function at  $Q^2=3.6$')
#
#  plt.show()

if __name__=='__main__':
    conf={}
    cwd = os.getcwd()
    #  conf['path2CJ'] ='../../external/CJLIB'
    conf['path2CT10'] ='../../external/PDF'
    conf['path2LSS']='../../external/LSSLIB'
    conf['path2DSS']='../../external/DSSLIB'
    
    conf['order']='LO'
    conf['aux']  =AUX()
    #  conf['_pdf'] =CJ(conf)
    conf['_pdf'] =CT10(conf)
    conf['_ppdf']=LSS(conf)
    conf['_ff']  =DSS(conf)
    
    conf['pdf']=PDF(conf)
    conf['ppdf']=PPDF(conf)
    conf['ff']=FF(conf)
    conf['gk']=GK(conf)
    
    stfuncs=STFUNCS(conf)
    x=0.25
    z=0.5
    Q2=2.4
    mu2=2.0
    E=11.0
    m=0.938
    y=Q2/(2*m*E*x)
    target='p'
    hadron='pi+' 
    
    bT = np.logspace(-3, 3, 60)
    FUUbT = [b*stfuncs.FUU_b(x,Q2,y,z,b,target,hadron) for b in bT]
    
    qT = np.linspace(0.01, Q2, 30)
    FUUqT = [stfuncs.FUU_q(x,Q2,y,z,q,target,hadron,0.3,0.5,Nmax = 20) for q in qT]
    FUUquad = [stfuncs.FUU_q_quad(x,Q2,y,z,q,target,hadron) for q in qT]
    
    ax = py.subplot(121)
    ax.plot(bT, FUUbT)
    ax.set_xlabel('b (1/GeV)', fontsize=10)
    ax.set_ylabel('FUU(b, x=0.25, z=0.5, Q2=2.4)', fontsize=10)
    ax.semilogx()
    
    ax = py.subplot(122)
    ax.errorbar(qT, FUUqT, [0]*len(qT),label = 'Ogata')
    ax.errorbar(qT, [FUUquad[i][0] for i in range(len(qT))], [FUUquad[i][1] for i in range(len(qT))], label = 'Quad')
    ax.set_xlabel('q (GeV)', fontsize=10)
    ax.set_ylabel('FUU(q, x=0.25, z=0.5, Q2=2.4)', fontsize=10)
    #ax.semilogx()
    py.tight_layout()
    ax.legend()
    py.show()