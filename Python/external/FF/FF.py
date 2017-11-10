#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad,dblquad,fixed_quad
#from line_profiler import LineProfiler
# DBC--removed factor multiplication 7/31/15

class FragFuncs(object):

  def __init__(self,fname):

    self.D={}
    self.load_table(fname)

  def load_table(self,fname):

    # read file
    F=open(fname,'r')
    L=F.readlines()
    F.close()

    # convert loaded file to floats
    L=[l.strip() for l in L]
    L=[[float(x) for x in l.split()] for l in L]

    # extract parts
    Q2=np.copy(L[0])
    X=np.copy(L[1])
    L=np.array(L[2:])

    # get number of partons in table
    npartons=len(L[0])

    # empy array for FF values
    FF=np.zeros((npartons,Q2.size,X.size))

    # fill array
    cnt=0
    for iX in range(X.size-1):
      for iQ2 in range(Q2.size):
        for iparton in range(npartons):
          if any([iparton==k for k in [0,1,2,6,7,8]]) : 
            factor=(1-X[iX])**4 * X[iX]**0.5
          elif any([iparton==k for k in [3,4]]) : 
            factor=(1-X[iX])**7 * X[iX]**0.3
          elif iparton==5: 
            factor=(1-X[iX])**4 * X[iX]**0.3
          FF[iparton,iQ2,iX]=L[cnt,iparton]#/factor
        cnt+=1

    LX=np.log(X)
    LQ2=np.log(Q2)

    D = self.D
    D['UTOT']=RectBivariateSpline(LQ2,LX,FF[0],kx=3,ky=3)
    D['DTOT']=RectBivariateSpline(LQ2,LX,FF[1],kx=3,ky=3)
    D['STOT']=RectBivariateSpline(LQ2,LX,FF[2],kx=3,ky=3)
    D['CTOT']=RectBivariateSpline(LQ2,LX,FF[3],kx=3,ky=3)
    D['BTOT']=RectBivariateSpline(LQ2,LX,FF[4],kx=3,ky=3)
    D['G']   =RectBivariateSpline(LQ2,LX,FF[5],kx=3,ky=3)
    D['UVAL']=RectBivariateSpline(LQ2,LX,FF[6],kx=3,ky=3)
    D['DVAL']=RectBivariateSpline(LQ2,LX,FF[7],kx=3,ky=3)
    D['SVAL']=RectBivariateSpline(LQ2,LX,FF[8],kx=3,ky=3)


    D['Up'] = lambda lQ2,lx: 0.5*(D['UTOT'](lQ2,lx)+D['UVAL'](lQ2,lx))
    D['UBp']= lambda lQ2,lx: 0.5*(D['UTOT'](lQ2,lx)-D['UVAL'](lQ2,lx))
    D['Dp'] = lambda lQ2,lx: 0.5*(D['DTOT'](lQ2,lx)+D['DVAL'](lQ2,lx))
    D['DBp']= lambda lQ2,lx: 0.5*(D['DTOT'](lQ2,lx)-D['DVAL'](lQ2,lx))
    D['Sp'] = lambda lQ2,lx: 0.5*(D['STOT'](lQ2,lx)+D['SVAL'](lQ2,lx))
    D['SBp']= lambda lQ2,lx: 0.5*(D['STOT'](lQ2,lx)-D['SVAL'](lQ2,lx))
    D['Cp'] = lambda lQ2,lx: 0.5*D['CTOT'](lQ2,lx)
    D['Bp'] = lambda lQ2,lx: 0.5*D['BTOT'](lQ2,lx)

  def get_xg(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    return self.D['G'](lQ2,lx)[0,0]#*(1-x)**4*x**0.3

  def get_xu(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['Up'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==-1:
      return D['UBp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5 *(D['Up'](lQ2,lx)[0,0]+D['UBp'](lQ2,lx)[0,0])#*(1-x)**4*x**0.5

  def get_xub(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['UBp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==-1:
      return D['Up'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Up'](lQ2,lx)[0,0]+D['UBp'](lQ2,lx)[0,0])#*(1-x)**4*x**0.5

  def get_xd(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['Dp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==-1:
      return D['DBp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Dp'](lQ2,lx)[0,0]+D['DBp'](lQ2,lx)[0,0])#*(1-x)**4*x**0.5

  def get_xdb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['DBp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==-1:
      return D['Dp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Dp'](lQ2,lx)[0,0]+D['DBp'](lQ2,lx)[0,0])#*(1-x)**4*x**0.5

  def get_xs(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['Sp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==-1:
      return D['SBp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Sp'](lQ2,lx)[0,0]+D['SBp'](lQ2,lx)[0,0])#*(1-x)**4*x**0.5

  def get_xsb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['SBp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==-1:
      return D['Sp'](lQ2,lx)[0,0]#*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Sp'](lQ2,lx)[0,0]+D['SBp'](lQ2,lx)[0,0])#*(1-x)**4*x**0.5

  def get_xc(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Cp'](lQ2,lx)[0,0]#*(1-x)**7*x**0.3

  def get_xcb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Cp'](lQ2,lx)[0,0]#*(1-x)**7*x**0.3

  def get_xb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Bp'](lQ2,lx)[0,0]#*(1-x)**7*x**0.3

  def get_xbb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Bp'](lQ2,lx)[0,0]#*(1-x)**7*x**0.3

  def get_FF(self,x,Q2,flav,charge):
    if   flav=='g':  return self.get_xg(x,Q2,charge)/x
    elif flav=='u':  return self.get_xu(x,Q2,charge)/x
    elif flav=='ub': return self.get_xub(x,Q2,charge)/x
    elif flav=='d':  return self.get_xd(x,Q2,charge)/x
    elif flav=='db': return self.get_xdb(x,Q2,charge)/x
    elif flav=='s':  return self.get_xs(x,Q2,charge)/x
    elif flav=='sb': return self.get_xsb(x,Q2,charge)/x
    elif flav=='c':  return self.get_xc(x,Q2,charge)/x
    elif flav=='cb': return self.get_xcb(x,Q2,charge)/x
    elif flav=='b':  return self.get_xb(x,Q2,charge)/x
    elif flav=='bb': return self.get_xbb(x,Q2,charge)/x


if __name__== "__main__":

  #import fDSS
  
  PI=FragFuncs('tables/PILO.TAB')
  ic=1
  
  #print PI.get_xg(0.6,10.0,ic)
  #print fDSS.fdss(1,ic,0,0.6,10.0)[8]
  zz = 0.05
  Q2 = 1.

  #zlist = []
  n1 = range(0,95)
  zlist = [zz + 0.01*nn for nn in n1]
  #vget = np.vectorize(self.get_FF)  
  FF = PI.get_FF
  #for zElem in zlist:
  #  u = FF(zElem,Q2,'u',ic)
  #  ub = FF(zElem,Q2,'ub',ic)
  #  d = FF(zElem,Q2,'d',ic)
  #  db = FF(zElem,Q2,'db',ic)
  #  s = FF(zElem,Q2,'s',ic)
  #  sb = FF(zElem,Q2,'sb',ic)
  #  c = FF(zElem,Q2,'c',ic)
  #  b = FF(zElem,Q2,'b',ic)
  #  g = FF(zElem,Q2,'g',ic)
  #  print '{0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(u,ub,d,db,s,sb,c,b,g)

  for zElem in zlist:
    u = PI.get_xu(zElem,Q2,ic)
    ub = PI.get_xub(zElem,Q2,ic)
    d = PI.get_xd(zElem,Q2,ic)
    db = PI.get_xdb(zElem,Q2,ic)
    s = PI.get_xs(zElem,Q2,ic)
    sb = PI.get_xsb(zElem,Q2,ic)
    c = PI.get_xc(zElem,Q2,ic)
    b = PI.get_xb(zElem,Q2,ic)
    g = PI.get_xg(zElem,Q2,ic)
    print '{0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(u,ub,d,db,s,sb,c,b,g)















