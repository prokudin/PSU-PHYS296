#!/usr/bin/env python
import sys,os
import numpy as np
from scipy.special import jv, jn_zeros, yv
import scipy.special as spec
from scipy.optimize import fsolve
from scipy.optimize import fmin
from scipy.integrate import quad
from scipy.optimize import fsolve

class OGATA:

    def __init__(self,xmin,xmax,nu):

        zero1 = jn_zeros(nu, 1)[0]
        h = fsolve(lambda h: xmin-zero1*np.tanh(np.pi/2*np.sinh(h/np.pi*zero1)), xmin)[0]
        k = fsolve(lambda k: xmax-np.pi*k*np.tanh(np.pi/2*np.sinh(h*k)), xmax)[0]
        N = int(k)
        #print 
        #print '\nogata N=',N
        #sys.exit()
        zeros=jn_zeros(nu,N)
        xi=zeros/np.pi
        Jp1=jv(nu+1,np.pi*xi)
        self.w=yv(nu, np.pi * xi) / Jp1
        get_psi=lambda t: t*np.tanh(np.pi/2*np.sinh(t))
        get_psip=lambda t:np.pi*t*(-np.tanh(np.pi*np.sinh(t)/2)**2 + 1)*np.cosh(t)/2 + np.tanh(np.pi*np.sinh(t)/2)
        self.knots=np.pi/h*get_psi(h*xi)
        self.Jnu=jv(nu,self.knots)
        self.psip=get_psip(h*xi)

    def invert(self,w, qT):
        F=w(self.knots/qT)/qT
        return 0.5*np.sum(self.w*F*self.Jnu*self.psip)


if __name__=='__main__':

  def Wtilde(bT,Q,sigma):
      M=1.0/Q
      V=1/sigma**2
      b=0.5*(-M+np.sqrt(M**2+4*V))
      a=V/b**2
      return bT**(a-1)*np.exp(-bT/b)/b**a/spec.gamma(a)
  
  def W(qT, Q, sigma, nu):
      M=1.0/Q
      V=1/sigma**2
      b=0.5*(-M+np.sqrt(M**2+4*V))
      a=V/b**2
      return 1/(2*np.pi)*spec.gamma(a+nu)/spec.gamma(a)*(b*qT/2.0)**nu*spec.hyp2f1((a+nu)/2.0, (a+nu+1.0)/2.0, nu+1.0, -qT**2.0*b**2.0)/spec.gamma(nu+1.0)

  nu=1
  Q=2.0
  sigma=1.0
  qTmin=0.1
  qTmax=3.0
  qT=np.linspace(qTmin,qTmax,10)
  xmin=1e-5*qTmin
  xmax=20*qTmax
  ogata=OGATA(xmin,xmax,1)

  for _qT in qT: 
    exact=W(_qT,Q,sigma,nu)
    approx=ogata.invert(lambda bT: Wtilde(bT,Q,sigma),_qT)
    rel=np.abs((exact-approx)/exact)*100
    print r'qT=%0.2f exact=%10.3e approx=%10.3e rel-err=%10.3f '%(_qT,exact,approx,rel)









