#!/usr/bin/env python
import sys,os
import numpy as np
import math
sys.path.append('./../../') #this appends/searches for the contents of the 
                  #the dir/folders that are two directories back
from tools.tools import load_config
from scipy.integrate import quad
from external.PDF.CT10 import CT10
#from external.CJLIB.CJ import CJ
from external.LSSLIB.LSS import LSS
from external.DSSLIB.DSS import DSS
from qcdlib.tmdlib import PDF,PPDF,FF,GK
from qcdlib.aux import AUX
import pylab as py
#import matplotlib.pyplot as plt
from tools.hankel.inverters import Ogata, Quad, Fix_Quad
import ogata

class STFUNCS:  # creating a class of

  def __init__(self,conf):
    self.aux=conf['aux']
    self.CF = conf['aux'].CF
    self.C2=1.0
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
    self.fquad = Fix_Quad()
    nu = 0
    qTmin=0.1
    qTmax=3.0
    xmin=0.5*qTmin
    xmax=10*qTmax
    self.ogata_fast=ogata.OGATA(xmin, xmax, nu)


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
    wq=self.get_wq(z,k1,k2,target,hadron) # +  self.get_wq_evolution(z,x,Q)
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

# bc
  def bc(self,b):
    return np.sqrt(b**2+(self.conf['gk'].bmin)**2)

# bstar
  def bstar(self,b):
    return b/np.sqrt(1+b**2/(self.conf['gk'].bmax)**2)

# mub
  def mub(self,b):
    return 2*np.exp(-np.euler_gamma)/self.bstar(b)

# gk function
#  def get_gk(self,b,z,x,Q):
#    wq=np.ones(len(self.e2))
#    return wq * self.conf['gk'].g2 * np.log(b/self.bstar(b)) * np.log( (z*Q)/(x*self.conf['gk'].Q0) )
#    return wq * self.conf['gk'].g2 * np.log(b/self.bstar(b)) * np.log( (Q)/(self.conf['gk'].Q0) )

    
    # gk function (b)
  def get_gk(self,b):
    wq=np.ones(len(self.e2))
    return wq * self.conf['gk'].g2 * np.log(b/self.bstar(b))

    
# Zeta prescription renormalization factor
  def ZETA_PRESCRIPTION(self, muf, zetaf, x, z, b):
    #CF = self.CF
    C0 = 2*np.exp(-np.euler_gamma)
    #gamma = self.conf['alphaS'].get_alphaS(muf**2)*CF/np.pi
    vf = 3.0/2.0+0.0
    #return gamma*np.log(muf*2*b/2./C0)*(np.log(zetaf*b/C0/muf)+vf)
    gk = self.get_gk(b)
    #pertutbative part of D
    #gamma0 = 4. * self.CF
    #beta0 = 11.-2./3.*3.  # self.nf
    #X=18.*self.conf['alphaS'].get_alphaS(muf**2)*np.log(muf*b*C0) 
    #print self.bstar(b)*C0
    #lX=np.log(1.-X)        
    #Dresummed=-8./27.*lX
    #return (gk-Dresummed)*(np.log(zetaf*b/C0/muf)+vf)
    return (gk)*(np.log(zetaf*b/C0/muf)+vf)
    
# intrinsic widths
  def get_width(self,b,z,k1,k2,target,hadron):
    return np.abs(self.conf[k1].widths[target])/4 + np.abs(self.conf[k2].widths[hadron])/(4*z**2)


# Structure functions in b space no evolution prescriprion
  def get_FX_b(self,i,x,z,Q2,pT,b,target,hadron):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    mu2=Q2
    if mu2>1000: mu2 = 1000.
    F=self.conf[k1].get_C(x,mu2,target)/(2*np.pi)
    D=self.conf[k2].get_C(z,mu2,hadron)/(2*np.pi*z**2)
    width=self.get_width(b,z,k1,k2,target,hadron)*b**2  
    K=self.get_K(i,x,Q2,z,pT,width,k1,k2,target,hadron)
    return 2*np.pi*np.sum(self.e2*K*F*D*np.exp(-width))  #sums up the contributions


# Structure functions in b space Zeta prescriprion
#  def get_FX_b(self,i,x,z,Q2,pT,b,target,hadron):
#    k1=self.D[i]['k1']
#    k2=self.D[i]['k2']
#    if k1==None or k2==None: return 0
#    Q = np.sqrt(Q2)
#    mu2=Q2
#    muf=Q
#    zetaf=Q2
#    if mu2>1000: mu2 = 1000.
#    F=self.conf[k1].get_C(x,mu2,target)/(2*np.pi)
#    D=self.conf[k2].get_C(z,mu2,hadron)/(2*np.pi*z**2)
#    width=self.get_width(b,z,k1,k2,target,hadron)*b**2  +  self.ZETA_PRESCRIPTION(muf, zetaf, x, z, b)
#    K=self.get_K(i,x,Q2,z,pT,width,k1,k2,target,hadron)
#    return 2*np.pi*np.sum(self.e2*K*F*D*np.exp(-width))  #sums up the contributions

# Perturbative evolution

  def gamma(self,mu):
    return conf['alphaS'].get_alphaS(mu**2)*self.CF/np.pi*(3.0/2.0 - 0.0)

  def gamma_k(self,mu):
    return 2*conf['alphaS'].get_alphaS(mu**2)*self.CF/np.pi

  def Int_gammas(self,mu1,mu2,Q):
    integrand=lambda mu: 1/mu*(2*self.gamma(mu)-2*np.log(Q/mu)*self.gamma_k(mu))
    return quad(integrand,mu1,mu2)[0]

  def Int_gammak(self,mu1,mu2):
    integrand=lambda mu: 1/mu*self.gamma_k(mu)
    return quad(integrand,mu1,mu2)[0]

  def Int_gammas_ope(self,bT,Q):
    return self.Int_gammas(self.mub(bT),self.C2*Q,Q) 

  def ope_evo(self,bT,Q):
    return self.Int_gammas_ope(bT,Q)

# Structure functions in b space CSS
  def get_FX_b_css(self,i,x,z,Q2,pT,b,target,hadron):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    #mu2=(self.mub(self.bc(b)))**2
    mu2=(self.mub(b))**2
    Q = np.sqrt(Q2)
    if mu2>1000: mu2 = 1000.
    #F=self.conf[k1].get_C(x,mu2,target)/(2*np.pi)
    #D=self.conf[k2].get_C(z,mu2,hadron)/(2*np.pi*z**2)
    F=self.conf[k1].get_ope_C(x,b,Q**2,Q,target,order=1)/(2*np.pi)
    D=self.conf[k2].get_ope_C(z,b,Q**2,Q,hadron,order=1)/(2*np.pi*z**2)
    width=self.get_width(b,z,k1,k2,target,hadron)*b**2\
          +self.get_gk(b)*np.log((Q)/(self.conf['gk'].Q0))\
          +self.ope_evo(b,Q)
    K=self.get_K(i,x,Q2,z,pT,width,k1,k2,target,hadron)
    return 2*np.pi*np.sum(self.e2*K*F*D*np.exp(-width))  #sums up the contributions
    
# Method by Jianwei
  def get_FX_b_jianwei(self,i,x,z,Q2,pT,b,target,hadron):
      if ( b > self.conf['gk'].b_max):
          return self.get_FX_b(i,x,z,Q2,pT,b,target,hadron)
      else:
          # Andrea and Tianbo will work here
          return 1./2.

# Structure function FUU in b space
  def FUU_b(self,x,Q2,y,z,q,b,target,hadron):
    factor = 1.0
    return factor*(self.get_FX_b(1,x,z,Q2,q,b,target,hadron))

  def FUU_q(self,x,Q2,y,z,q,target,hadron,Nmax = 13):
    nu = 0
    Q = np.sqrt(Q2)
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron))
    return 2*np.pi*self.ogata.adog3(w, q, nu, Nmax, Q)

  def FUU_q_quad(self,x,Q2,y,z,q,target,hadron,eps = 1e-3):
    nu = 0
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron))
    inv = self.quad.quadinv(w, q, nu, eps)
    return 2*np.pi*inv[0], 2*np.pi*inv[1]

  def FUU_q_fquad(self,x,Q2,y,z,q,target,hadron,num):
    nu = 0
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron))
    inv = self.fquad.fix_quadinv(w, q, nu, num)
    return 2*np.pi*inv#[0], 2*np.pi*inv[1]

  def FUU_fast(self,x,Q2,y,z,q,target,hadron):
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron))
    return 2*np.pi*self.ogata_fast.invert(w,q)


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
    import qcdlib.aux
    import qcdlib.alphaS
    
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
    
    conf['gk']=GK(conf)
    conf['pdf']=PDF(conf)
    conf['ppdf']=PPDF(conf)
    conf['ff']=FF(conf)
    conf['Q20'] = conf['gk'].Q0**2
    conf['alphaSmode']='backward'
    conf['alphaS']=qcdlib.alphaS.ALPHAS(conf)
    
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
    
    bT = np.logspace(-1, 0, 30)
    pT=1.0
    FUUbCSS = [b*stfuncs.get_FX_b_css(1,x,z,Q2,pT,b,target,hadron) for b in bT]
    FUUb = [b*stfuncs.get_FX_b(1,x,z,Q2,pT,b,target,hadron) for b in bT]
    ratio = [FUUbCSS[i]/FUUb[i] for i in range(len(bT))]
    
    ax = py.subplot(121)
    ax.plot(bT, FUUbCSS, label = 'CSS')
    ax.plot(bT, FUUb, label = 'Gauss')
    ax.set_xlabel('b_T', fontsize=10)
    ax.set_ylabel('FUU(b, x='+str(x)+', z='+str(z)+', Q2='+str(Q2)+')', fontsize=10)
    ax.semilogx()
    #ax.semilogy()
    ax.legend()
    ax=py.subplot(122)
    ax.plot(bT,ratio,label='ratio')
    
    ax.set_xlabel('b_T', fontsize=10)
    ax.set_ylabel('CSS/Gauss(b, x='+str(x)+', z='+str(z)+', Q2='+str(Q2)+')', fontsize=10)
    #ax.set_ylabel('pert', fontsize=10)
    ax.set_ylim(0, 2)
    ax.semilogx()
    
    py.tight_layout()
    ax.legend()
    py.show()
