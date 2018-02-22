#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../') 
import numpy as np
import time
from tools.tools import load_config
from scipy.integrate import quad
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

  def get_params(self):
    pre_params=load_config(self.conf['path2upol_compass'])
    params = {'pdf': ['widths0 valence', 'widths0 sea'], 'ff': ['widths0 pi+ unfav', 'widths0 k+ fav', 'widths0 pi+ fav', 'widths0 k+ unfav'], 'gk': ['Q0', 'bmax', 'g2', 'bmin']}
    newparams = {'pdf': {}, 'ff': {}, 'gk': {}}
    for key in params.keys():
      for param in params[str(key)]:
        newparams[str(key)][str(param)] = pre_params['params'][str(key)][str(param)]['value']
    return newparams['pdf']


  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
    self.setup()

  def set_default_params(self):

    self.widths0={}
    self.get_params()['widths0 valence']
    self.widths0['valence']=self.get_params()['widths0 valence'] #0.3
    self.widths0['sea']= self.get_params()['widths0 sea'] #0.3

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

class FF(CORE):

  def get_params(self):
    pre_params=load_config(self.conf['path2upol_compass'])
    params = {'pdf': ['widths0 valence', 'widths0 sea'], 'ff': ['widths0 pi+ unfav', 'widths0 k+ fav', 'widths0 pi+ fav', 'widths0 k+ unfav'], 'gk': ['Q0', 'bmax', 'g2', 'bmin']}
    newparams = {'pdf': {}, 'ff': {}, 'gk': {}}
    for key in params.keys():
      for param in params[str(key)]:
        newparams[str(key)][str(param)] = pre_params['params'][str(key)][str(param)]['value']
    return newparams['ff']

  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
    self.setup()

  def set_default_params(self):

    self.widths0={}
    self.widths0['pi+ fav']  =self.get_params()['widths0 pi+ fav']   #0.12
    self.widths0['pi+ unfav']=self.get_params()['widths0 pi+ unfav'] #0.12
    self.widths0['k+ fav']   =self.get_params()['widths0 k+ fav']    #0.12
    self.widths0['k+ unfav'] =self.get_params()['widths0 k+ unfav']  #0.12
    
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

# class for "toy" evolution
class GK(CORE):

  def get_params(self):
    pre_params=load_config(self.conf['path2upol_compass'])
    params = {'pdf': ['widths0 valence', 'widths0 sea'], 'ff': ['widths0 pi+ unfav', 'widths0 k+ fav', 'widths0 pi+ fav', 'widths0 k+ unfav'], 'gk': ['Q0', 'bmax', 'g2', 'bmin']}
    newparams = {'pdf': {}, 'ff': {}, 'gk': {}}
    for key in params.keys():
      for param in params[str(key)]:
        newparams[str(key)][str(param)] = pre_params['params'][str(key)][str(param)]['value']
    return newparams['gk']

  def __init__(self,conf):
    self.aux=conf['aux']
    self.conf=conf
    self.set_default_params()
#    self.setup()

  def set_default_params(self):
    self.g2=self.get_params()['g2']#0.1
    self.Q0=self.get_params()['Q0']#1.0
    self.bmax=self.get_params()['bmax']#1.0
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

  conf={}
#  conf['path2CJ'] ='../external/CJLIB'
  conf['path2upol_compass'] ='../fitlab/inputs/upol_compass.py'
  conf['path2CT10'] ='../external/PDF'
  conf['path2LSS']='../external/LSSLIB'
  conf['path2GRSV']='../external/GRSVLIB'
  conf['path2DSS']='../external/DSSLIB'

  conf['order']='LO'
  conf['aux']=AUX()
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
  dist=['pdf','ppdf','ff']
  for k in dist:
    print k
    print conf[k].get_C(x,Q2)





