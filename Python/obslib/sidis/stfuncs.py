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
from tools.config import conf

class STFUNCS:  # creating a class of

  def __init__(self):
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
    return z**2*np.abs(conf[k1].widths[target]) + np.abs(conf[k2].widths[hadron])

  def get_wq_evolution(self,z,x,Q):
    wq=np.ones(len(self.e2))
    return wq * 4*z**2 * conf['gk'].g2 * np.log( (z*Q)/(x*conf['gk'].Q0) )
#    return wq * 4*z**2 * conf['gk'].g2 * np.log( (Q)/(conf['gk'].Q0) )

  def get_gauss(self,z,pT,wq):
    return np.exp(-pT**2/wq)/(np.pi*wq)

# Structure functions in pT space
  def get_FX(self,i,x,z,Q2,pT,target,hadron):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    mu2=Q2
    Q = np.sqrt(Q2)
    F=conf[k1].get_C(x,mu2,target)
    D=conf[k2].get_C(z,mu2,hadron)
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
    return np.sqrt(b**2+(conf['gk'].bmin)**2)

# bstar
  def bstar(self,b):
    return b/np.sqrt(1+b**2/(conf['gk'].bmax)**2)

# mub
  def mub(self,b):
    return 2*np.exp(-np.euler_gamma)/self.bstar(b)

# gk function
#  def get_gk(self,b,z,x,Q):
#    wq=np.ones(len(self.e2))
#    return wq * conf['gk'].g2 * np.log(b/self.bstar(b)) * np.log( (z*Q)/(x*conf['gk'].Q0) )
#    return wq * conf['gk'].g2 * np.log(b/self.bstar(b)) * np.log( (Q)/(conf['gk'].Q0) )

    
    # gk function (b)
  def get_gk(self,b):
    wq=np.ones(len(self.e2))
    return wq * conf['gk'].g2 * np.log(b/self.bstar(b))

    
# Zeta prescription renormalization factor
  def ZETA_PRESCRIPTION(self, muf, zetaf, x, z, b):
    #CF = self.CF
    C0 = 2*np.exp(-np.euler_gamma)
    #gamma = conf['alphaS'].get_alphaS(muf**2)*CF/np.pi
    vf = 3.0/2.0+0.0
    #return gamma*np.log(muf*2*b/2./C0)*(np.log(zetaf*b/C0/muf)+vf)
    gk = self.get_gk(b)
    #pertutbative part of D
    #gamma0 = 4. * self.CF
    #beta0 = 11.-2./3.*3.  # self.nf
    #X=18.*conf['alphaS'].get_alphaS(muf**2)*np.log(muf*b*C0) 
    #print self.bstar(b)*C0
    #lX=np.log(1.-X)        
    #Dresummed=-8./27.*lX
    #return (gk-Dresummed)*(np.log(zetaf*b/C0/muf)+vf)
    return (gk)*(np.log(zetaf*b/C0/muf)+vf)
    
# intrinsic widths
  def get_width(self,b,z,k1,k2,target,hadron):
    return np.abs(conf[k1].widths[target])/4 + np.abs(conf[k2].widths[hadron])/(4*z**2)


# Structure functions in b space no evolution prescriprion
  def get_FX_b(self,i,x,z,Q2,pT,b,target,hadron,order):
    k1=self.D[i]['k1']
    k2=self.D[i]['k2']
    if k1==None or k2==None: return 0
    mu2=Q2
    if mu2>1000: mu2 = 1000.
    F=conf[k1].get_C(x,mu2,target)/(2*np.pi)
    D=conf[k2].get_C(z,mu2,hadron)/(2*np.pi*z**2)
    width=self.get_width(b,z,k1,k2,target,hadron)*b**2  
    K=self.get_K(i,x,Q2,z,pT,width,k1,k2,target,hadron)
    return 2*np.pi*np.sum(self.e2*K*F*D*np.exp(-width))  #sums up the contributions


# Structure functions in b space Zeta prescriprion
#  def get_FX_b(self,i,x,z,Q2,pT,b,target,hadron,order):
#    k1=self.D[i]['k1']
#    k2=self.D[i]['k2']
#    if k1==None or k2==None: return 0
#    Q = np.sqrt(Q2)
#    mu2=Q2
#    muf=Q
#    zetaf=Q2
#    if mu2>1000: mu2 = 1000.
#    F=conf[k1].get_C(x,mu2,target)/(2*np.pi)
#    D=conf[k2].get_C(z,mu2,hadron)/(2*np.pi*z**2)
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

  def ope_evo(self,bT,Q,order):
    if order == 0: return 0
    elif order == 1: return self.Int_gammas_ope(bT,Q)-2*self.Int_gammak(self.mub(bT),self.C2*conf['gk'].Q0)*np.log(Q/conf['gk'].Q0)

# Hard Factor

  def get_H(self,Q,order):
    if order == 0: return 1
    elif order == 1:
      alphaS=conf['alphaS'].get_alphaS(Q**2)
    #return 1+4*self.CF*(-3*np.log(self.C2)-0.5*np.log(self.C2**2)**2-4)*alphaS/(4*np.pi)
      return 1+4*self.CF*(3*np.log(self.C2)-0.5*np.log(self.C2**2)**2-4)*alphaS/(4*np.pi)

# Structure functions in b space CSS
#  def get_FX_b_css(self,i,x,z,Q2,pT,b,target,hadron,order):
#    k1=self.D[i]['k1']
#    k2=self.D[i]['k2']
#    if k1==None or k2==None: return 0
#    #mu2=(self.mub(self.bc(b)))**2
#    mu2=(self.mub(b))**2
#    Q = np.sqrt(Q2)
#    if mu2>1000: mu2 = 1000.
#    #F=conf[k1].get_C(x,mu2,target)/(2*np.pi)
#    #D=conf[k2].get_C(z,mu2,hadron)/(2*np.pi*z**2)
#    F=conf[k1].get_ope_C(x,self.bstar(b),mu2,np.sqrt(mu2),target,order)/(2*np.pi)
#    D=conf[k2].get_ope_C(z,self.bstar(b),mu2,np.sqrt(mu2),hadron,order)/(2*np.pi*z**2)
#    width=self.get_width(b,z,k1,k2,target,hadron)*b**2\
#          +self.get_gk(b)*np.log((Q)/(conf['gk'].Q0))\
#          -self.ope_evo(b,Q,order)
#    K=self.get_K(i,x,Q2,z,pT,width,k1,k2,target,hadron)
#    return 2*np.pi*np.sum(self.e2*K*F*D*np.exp(-width))  #sums up the contributions
    
# Method by Jianwei
  def get_FX_b_jianwei(self,i,x,z,Q2,pT,b,target,hadron):
      if ( b > conf['gk'].b_max):
          return self.get_FX_b(i,x,z,Q2,pT,b,target,hadron)
      else:
          # Andrea and Tianbo will work here
          return 1./2.

# Structure function FUU in b space
  def FUU_b(self,x,Q2,y,z,q,b,target,hadron,order=0):
    Q = np.sqrt(Q2)
    H = self.get_H(Q,order)
    return H*(self.get_FX_b(1,x,z,Q2,q,b,target,hadron,order))

  def FUU_q(self,x,Q2,y,z,q,target,hadron,order,Nmax = 13):
    nu = 0
    Q = np.sqrt(Q2)
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron,order))
    return 2*np.pi*self.ogata.adog3(w, q, nu, Nmax, Q)

  def FUU_q_quad(self,x,Q2,y,z,q,target,hadron,order,eps = 1e-3):
    nu = 0
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron,order))
    inv = self.quad.quadinv(w, q, nu, eps)
    return 2*np.pi*inv[0], 2*np.pi*inv[1]

  def FUU_q_fquad(self,x,Q2,y,z,q,target,hadron,num,order):
    nu = 0
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron,order))
    inv = self.fquad.fix_quadinv(w, q, nu, num)
    return 2*np.pi*inv#[0], 2*np.pi*inv[1]

  def FUU_fast(self,x,Q2,y,z,q,target,hadron,order):
    w = np.vectorize(lambda b: b*self.FUU_b(x,Q2,y,z,q,b,target,hadron,order))
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
    
    # Compare CSS and Gaussian in b space
    bT = np.logspace(-2, 0.7, 30)
    pT=1.0
    order=0
    ordercss=1
    FUUbCSS = [b*stfuncs.get_FX_b_css(1,x,z,Q2,pT,b,target,hadron,ordercss) for b in bT]
    FUUb = [b*stfuncs.get_FX_b(1,x,z,Q2,pT,b,target,hadron,order) for b in bT]
    evo = [stfuncs.ope_evo(b,np.sqrt(Q2),ordercss) for b in bT]
    ratio = [FUUbCSS[i]/FUUb[i] for i in range(len(bT))]
    
    ax = py.subplot(121)
    ax.plot(bT, FUUbCSS, label = 'CSS')
    ax.plot(bT, FUUb, label = 'Gauss')
    ax.set_xlabel('b', fontsize=10)
    ax.set_ylabel('b FUU(b, x='+str(x)+', z='+str(z)+', Q2='+str(Q2)+')', fontsize=10)
    ax.semilogx()
    ax.semilogy()
    ax.legend()
#    ax=py.subplot(132)
#    ax.plot(bT,ratio,label='ratio')
#    
#    ax.set_xlabel('b_T', fontsize=10)
#    ax.set_ylabel('CSS/Gauss(b, x='+str(x)+', z='+str(z)+', Q2='+str(Q2)+')', fontsize=10)
#    #ax.set_ylabel('pert', fontsize=10)
#    ax.set_ylim(0, 2)
#    ax.semilogx()
#    
#    ax = py.subplot(133)
#    ax.plot(bT,evo, label='pert')
#    ax.set_xlabel('bT', fontsize=10)
#    ax.set_ylabel('Sudakov', fontsize = 10)
#    ax.semilogx()
#    
#    py.tight_layout()
#    ax.legend()
#    py.show()

# Compare gaussian pdf, ff to css in bspace 
    k1=stfuncs.D[1]['k1']
    k2=stfuncs.D[1]['k2']
    bT = np.logspace(-2, 1, 30)
    Q = np.sqrt(Q2)
    zeta = lambda b: stfuncs.mub(b)**2
    mu = lambda b: stfuncs.mub(b)
    opeff = [stfuncs.conf[k1].get_ope_C(x,stfuncs.bstar(b),zeta(b),mu(b),target,1)[1] for b in bT]
    gaussff = [stfuncs.conf[k1].get_C(x,mu2,target)[1] for b in bT]
    opepdf = [stfuncs.conf[k2].get_ope_C(x,stfuncs.bstar(b),zeta(b),mu(b),hadron,1)[1] for b in bT]
    gausspdf = [stfuncs.conf[k2].get_C(x,mu2,hadron)[1] for b in bT]
    ffratio = [opeff[i]/gaussff[i] for i in range(len(bT))]
    pdfratio = [opepdf[i]/gausspdf[i] for i in range(len(bT))]
    prod = [ffratio[i]*pdfratio[i] for i in range(len(bT))]
    
    ax = py.subplot(122)
    ax.plot(bT,ffratio,label = 'FF')
    ax.plot(bT,pdfratio,label='PDF')
    ax.plot(bT,prod,label='FF*PDF')
    ax.set_xlabel('b')
    ax.set_ylabel(r'$TMD_{CSS}/TMD_{Gauss}$(b,x=0.25,z=0.5,Q2=2.4)')
    ax.semilogx()
    ax.legend()
    #ax.set_ylim(-2, 2)
    py.tight_layout()
    py.show()
