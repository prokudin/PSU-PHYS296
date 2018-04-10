#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../') 
import numpy as np
import time
from scipy.integrate import quad,fixed_quad
#from external.CJLIB.CJ import CJ
from external.PDF.CT10 import CT10
from external.LSSLIB.LSS import LSS
from external.DSSLIB.DSS import DSS
from external.GRSVLIB.GRSV import GRSV
from aux import AUX
from scipy.special import gamma

class CORE:

  def p2n(self,p):
    n=np.copy(p)
    u  =n[1]
    ub =n[2]
    d  =n[3]
    db =n[4]
    n[1] = d
    n[2] = db
    n[3] = u
    n[4] = ub
    return n

  def pip2pim(self,pip):
    pim=np.copy(pip)
    u  =pip[1]
    ub =pip[2]
    d  =pip[3]
    db =pip[4]
    pim[1] = ub
    pim[2] = u
    pim[3] = db
    pim[4] = d
    return pim

  def kp2km(self,kp):
    km=np.copy(kp)
    u  =kp[1]
    ub =kp[2]
    s  =kp[5]
    sb =kp[6]
    km[1] = ub
    km[2] = u
    km[5] = sb
    km[6] = s
    return km

  def beta(self,a,b):
    return gamma(a)*gamma(b)/gamma(a+b)

  def get_shape(self,x,p):
    if self.conf['shape']==0:
       return  p[0]*x**p[1]*(1-x)**p[2]*(1+p[3]*x+p[4]*x**2)
    elif self.conf['shape']==1:
       norm=self.beta(1+p[1],p[2]+1)+p[3]*self.beta(1+p[1]+1,p[2]+1)+p[4]*self.beta(1+p[1]+2,p[2]+1)
       return  p[0]*x**p[1]*(1-x)**p[2]*(1+p[3]*x+p[4]*x**2)/norm

  def get_collinear(self,x,hadron):
    N=np.zeros(11)
    for i in range(11): 
      N[i]=self.get_shape(x,self.shape[hadron][i])
    return N

  def get_gauss(self,kT2,hadron):
    return np.exp(-kT2/self.widths[hadron])/np.pi/self.widths[hadron] 

  def get_tmd(self,x,Q2,kT2,hadron):
    C=self.get_C(x,Q2,hadron)
    gauss=self.get_gauss(kT2,hadron)
    return self.K[hadron]*C*gauss

class PDF(CORE):
    
  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
    self.setup()
    self.forder=[0,3,1,5,7,9,10,8,6,2,4]
    self.method='quad'
    self.CF = conf['aux'].CF
    self.TF=conf['aux'].TF
    self.gamma_E=conf['aux'].euler
    self.C1=2*np.exp(-self.gamma_E)
    self.bmax=self.conf['gk'].bmax

  # aux functions
  def integrator(self,f,xmin,xmax):
    if self.method=='quad':
      return quad(f,xmin,xmax)[0]
    elif self.method=='gauss':
      return fixed_quad(np.vectorize(f),xmin,xmax,n=200)[0]

  def plus(self,x,xh,f):
    return  2/(1-xh)*(f(xh)-f(1)) + 2*np.log(1-x)*f(1)/(1-x)

  def bstar(self,bT,Q):
    return bT/np.sqrt(1+bT**2/self.bmax**2)

  def mub(self,bT,Q):
    return self.C1/self.bstar(bT,Q)

  def set_default_params(self):

    self.widths0={}
    self.widths0['valence']=0.3
    self.widths0['sea']=0.3

    self.widths={}
    self.widths['p']=np.ones(11)

    self.K={}
    self.K['p']=np.ones(11)
    self.K['n']=np.ones(11)

  def setup(self):
    for i in range(11):
      if i==1 or i==3:
        self.widths['p'][i]=self.widths0['valence']
      else:
        self.widths['p'][i]=self.widths0['sea']

    self.widths['n']=self.p2n(self.widths['p'])

  def get_C(self,x,Q2,target='p'):
    C=self.conf['_pdf'].get_f(x,Q2)
    C[0]=0 # glue is not supported
    if target=='n': C=self.p2n(C)
    return C

  def get_C_ope(self,i,x,Q2,target='p'):
    C=self.conf['_pdf'].get_f(x,Q2)#[self.forder]
    if target=='n': C=self.p2n(C)
    return C[i]

  def integrand_jpj(self,x,xh,f,bT,zetaF,mu,order):
    integrand=f(x)/(1-x)
    if order>0:
      alphaS=conf['alphaS'].get_alphaS(mu**2)
      expression = 2*(np.log(2.0/(mu*bT))-self.gamma_E)*(self.plus(x,xh,lambda xh:f(x/xh)/xh)\
            -(1+xh)*f(x/xh)/xh)+ (1-xh)/xh*f(x/xh)\
          +( -0.5*(np.log(bT**2*mu**2)-2*(np.log(2)-self.gamma_E))**2\
              -(np.log(bT**2*mu**2) - 2*(np.log(2)-self.gamma_E)) * np.log(zetaF/mu**2))*f(x)/(1-x)
      integrand+=alphaS*self.CF/2.0/np.pi*expression
    return integrand

  def integrand_jg(self,x,xh,f,bT,zetaF,mu,order):
    integrand=0
    if order>0:
      alphaS=conf['alphaS'].get_alphaS(mu**2)
      integrand+=alphaS*self.TF/2.0/np.pi*(2*(1-2*xh*(1-xh))*(np.log(2.0/bT/mu)\
        -self.gamma_E)+2*xh*(1-xh))/xh*f(x/xh)
    return integrand

  def get_ope_C(self,x,bT,zetaF,mu,target='p',order=0):
    opelst = [0]*len(self.forder)
    for i in range(len(self.forder)):
      f=lambda xx: self.get_C_ope(i,xx,mu**2,target)
      if i!=0:   integrand=lambda xh: self.integrand_jpj(x,xh,f,bT,zetaF,mu,order)
      elif i==0: integrand=lambda xh: self.integrand_jg( x,xh,f,bT,zetaF,mu,order)
      opelst[i]=self.integrator(integrand,x,1)
    return opelst

class FF(CORE):
    
  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
    self.setup()
    self.forder=[0,3,1,5,7,9,10,8,6,2,4]
    self.method='quad'
    self.CF = conf['aux'].CF
    self.TF=conf['aux'].TF
    self.gamma_E=conf['aux'].euler
    self.C1=2*np.exp(-self.gamma_E)
    self.bmax=self.conf['gk'].bmax

  # aux functions
  def integrator(self,f,xmin,xmax):
    if self.method=='quad':
      return quad(f,xmin,xmax)[0]
    elif self.method=='gauss':
      return fixed_quad(np.vectorize(f),xmin,xmax,n=200)[0]

  def plus(self,x,xh,f):
    return  2/(1-xh)*(f(xh)-f(1)) + 2*np.log(1-x)*f(1)/(1-x)

  def bstar(self,bT,Q):
    return bT/np.sqrt(1+bT**2/self.bmax**2)

  def mub(self,bT,Q):
    return self.C1/self.bstar(bT,Q)

  def set_default_params(self):

    self.widths0={}
    self.widths0['pi+ fav']  =0.12
    self.widths0['pi+ unfav']=0.12
    self.widths0['k+ fav']   =0.12
    self.widths0['k+ unfav'] =0.12
    
    self.widths={}
    self.widths['pi+']=np.ones(11)
    self.widths['k+'] =np.ones(11)

    self.K={}
    self.K['pi+']=np.ones(11)
    self.K['k+'] =np.ones(11)
    self.K['pi-']=np.ones(11)
    self.K['k-'] =np.ones(11)

  def setup(self):
    # u  1
    # ub 2
    # d  3
    # db 4
    # s  5
    # sb 6
    # c  7
    # cb 8
    # b  9
    # bb 10
    for i in range(1,11):
      if i==1 or i==4:
        self.widths['pi+'][i]=np.copy(self.widths0['pi+ fav'])
      else: 
        self.widths['pi+'][i]=np.copy(self.widths0['pi+ unfav'])

    for i in range(1,11):
      if i==1 or i==6:
        self.widths['k+'][i]=np.copy(self.widths0['k+ fav'])
      else: 
        self.widths['k+'][i]=np.copy(self.widths0['k+ unfav'])

    self.widths['pi-']=self.pip2pim(self.widths['pi+'])
    self.widths['k-'] =self.kp2km(self.widths['k+'])

  def get_C(self,x,Q2,hadron='pi+'):
    C=self.conf['_ff'].get_f(x,Q2,hadron)
    C[0]=0 # glue is not supported
    return C

  def get_C_ope(self,i,x,Q2,hadron='pi+'):
    C=self.conf['_ff'].get_f(x,Q2,hadron)
    return C[i]

  def integrand_jjp(self,z,zh,f,bT,zetaD,mu,order):
    integrand=1/(1-z)*f(z)/z**2
    if order>0:
      alphaS=conf['alphaS'].get_alphaS(mu**2)
      expression = self.plus(z,zh,lambda zh:2*(np.log(2.0*zh/(mu*bT))-self.gamma_E)*f(z/zh)*zh)\
          +2*(np.log(2.0*zh/(mu*bT))-self.gamma_E)*(1/zh**2+1/zh)*zh*f(z/zh)\
          +(1/zh**2 - 1/zh)*zh*f(z/zh) + 1/(1-z)*(-0.5*(np.log(bT**2*mu**2)-2*(np.log(2)-self.gamma_E))**2\
          -(np.log(bT**2*mu**2) - 2*(np.log(2)-self.gamma_E)) * np.log(zetaD/mu**2))*f(z)
      integrand+=alphaS*self.CF/2.0/np.pi*expression/z**2
    return integrand

  def integrand_gj(self,z,zh,f,bT,zetaD,mu,order):
    integrand=0
    if order>0:
      alphaS=conf['alphaS'].get_alphaS(mu**2)
      integrand+=alphaS*self.CF/2.0/np.pi/zh**3*(2*(1+(1-zh)**2)*(np.log(2.0*zh/bT/mu)\
        -self.gamma_E)+zh**2)*zh*f(z/zh)/z**2
    return integrand
  
  def get_ope_C(self,z,bT,zetaD,mu,hadron='pi+',order=0):
    opelst = [0]*len(self.forder)
    for i in range(len(self.forder)):
        f=lambda zz: self.get_C_ope(i,zz,mu**2,hadron)
        if i!=0:   integrand=lambda zh:self.integrand_jjp(z,zh,f,bT,zetaD,mu,order)
        elif i==0: integrand=lambda zh:self.integrand_gj( z,zh,f,bT,zetaD,mu,order)
        opelst[i]=z**2*self.integrator(integrand,z,1)
    return opelst

# class for "toy" evolution
class GK(CORE):

  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
#    self.setup()

  def set_default_params(self):
    self.g2=0.1
    self.Q0=1.
    self.bmax=1.
    self.bmin=0.1

#  def setup(self):


class PPDF(CORE):

  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
    self.setup()

  def set_default_params(self):
    self.widths0={}
    self.widths0['valence']=0.19
    self.widths0['sea']=0.19

    self.widths={}
    self.widths['p']=np.ones(11)

    self.K={}
    self.K['p']=np.ones(11)
    self.K['n']=np.ones(11)

  def setup(self):
    for i in range(11):
      if i==1 or i==3:
        self.widths['p'][i]=self.widths0['valence']
      else:
        self.widths['p'][i]=self.widths0['sea']

    self.widths['n']=self.p2n(self.widths['p'])

  def get_C(self,x,Q2,target='p'):
    C=self.conf['_ppdf'].get_f(x,Q2)
    C[0]=0 # glue is not supported
    if target=='n': C=self.p2n(C)
    return C



if __name__=='__main__':

  import qcdlib.alphaS
  conf={}
#  conf['path2CJ'] ='../external/CJLIB'
  conf['path2CT10'] ='../external/PDF'
  conf['path2LSS']='../external/LSSLIB'
  conf['path2GRSV']='../external/GRSVLIB'
  conf['path2DSS']='../external/DSSLIB'

  conf['order']='LO'
  conf['Q20']=1.0
  conf['alphaSmode']='backward'
  conf['aux']=AUX()
  conf['alphaS']=qcdlib.alphaS.ALPHAS(conf)
  conf['gk']=GK(conf)
#  conf['_pdf']=CJ(conf)
  conf['_pdf']=CT10(conf)
#  conf['_ppdf']=LSS(conf)
  conf['_ppdf']=GRSV(conf)
  conf['_ff']=DSS(conf)
  conf['hadron']='pi'

  conf['pdf']=PDF(conf)
  conf['ppdf']=PPDF(conf)
  conf['ff']=FF(conf)
  conf['gk']=GK(conf)

  x=0.15
  Q2=2.4
  dist=['pdf','ff']
  for k in dist:
    print k
    print conf[k].get_C(x,Q2)
    b=1.0
    Q = np.sqrt(Q2)
    mu=conf[k].mub(b,np.sqrt(Q))
    zeta = mu**2
    print conf[k].get_ope_C(x,b,zeta,mu,order=1)





