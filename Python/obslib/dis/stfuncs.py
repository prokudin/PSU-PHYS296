#!/usr/bin/env python
import numpy as np
import pylab as py
from scipy.special import gamma
from scipy.integrate import quad as _quad
from scipy.integrate import fixed_quad
from qcdlib.aux import AUX
from qcdlib.alphaS import ALPHAS
#from external.CJLIB.CJ import CJ
from external.PDF.CT10 import CT10

class STFUNCS:
  
  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.CF = self.aux.CF
    self.TR = self.aux.TR
    self.CA = self.aux.CA
    self.eU2=4.0/9.0
    self.eD2=1.0/9.0
    self.storage={}
    
  def integrator(self,f,xmin,xmax,method,n=200):
    f=np.vectorize(f)
    if method=='quad':
      return quad(f,xmin,xmax)[0]
    elif method=='gauss':
      return fixed_quad(f,xmin,xmax,n=n)[0]
    
  def log_plus(self,z,f,x):
    return np.log(1-z)/(1-z)*(f(x/z)/z-f(x)) + 0.5*np.log(1-x)**2*f(x)/(1-x)

  def one_plus(self,z,f,x):
    return 1/(1-z)*(f(x/z)/z-f(x))+ np.log(1-x)*f(x)/(1-x)
    
  def Cq(self,z,f,x):
    return self.CF*(\
        2*self.log_plus(z,f,x)-1.5*self.one_plus(z,f,x)\
      +(-(1+z)*np.log(1-z)-(1+z*z)/(1-z)*np.log(z)+3+2*z)*f(x/z)/z\
      -(np.pi**2/3+4.5)*f(x)/(1-x))
    
  def Cg(self,z,f,x):
    return self.TR*(((1-z)**2+z*z)*np.log((1-z)/z)-8*z*z+8*z-1)*f(x/z)    
 
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
 
  def qplus(self,PDF):
    if self.hadron=='n': PDF=self.p2n(PDF)
    out=self.eU2*(PDF[1]+PDF[2])\
       +self.eD2*(PDF[3]+PDF[4])\
       +self.eD2*(PDF[5]+PDF[6])
    if self.Nf>3:
      out+=self.eU2*(PDF[7]+PDF[8])
    if self.Nf>4:
      out+=self.eD2*(PDF[9]+PDF[10])
    return out
      
  def glue(self,PDF):
    factor=self.eU2 + 2*self.eD2
    if self.Nf>3:
      factor+=self.eD2
    if self.Nf>4:
      factor+=self.eU2
    return 2*factor*PDF[0]
      
  def integrand(self,x,z,Q2):
    return self.Cq(z,lambda y:self.qplus(self.get_PDFs(y)),x)\
         + self.Cg(z,lambda y:self.glue(self.get_PDFs(y)),x)

  def get_F2(self,x,Q2,hadron,method='gauss',n=20):
    if (x,Q2,hadron) not in self.storage:
      self.hadron=hadron
      self.Nf=self.conf['alphaS'].get_Nf(Q2)
      #alpi=self.conf['alphaS'].get_alphaS(Q2)/(2.*np.pi)
      self.get_PDFs=lambda y: self.conf['pdf-NLO'].get_f(y,Q2)
      LO=self.qplus(self.get_PDFs(x))
      #integrand=lambda z:self.integrand(x,z,Q2)
      #NLO=self.integrator(integrand,x,1,method,n=n)
      self.storage[(x,Q2,hadron)]=x*(LO)#+alpi*NLO)
    return self.storage[(x,Q2,hadron)]

if __name__=='__main__':

  conf={}
  conf['alphaSmode']='backward'
  conf['order']='LO'
  conf['Q20']=1
  conf['aux']=AUX()
  conf['alphaS']=ALPHAS(conf)
#  conf['path2CJ'] ='../../external/CJLIB'
#  conf['pdf-NLO']=CJ(conf)
  conf['path2CT10'] ='../../external/PDF'
  conf['pdf-NLO']=CT10(conf)
  stfuncs=STFUNCS(conf)

  x=0.5
  Q2=1000.
  print stfuncs.get_F2(x,Q2,'p',method='gauss')
  print stfuncs.get_F2(x,Q2,'n',method='gauss')

























