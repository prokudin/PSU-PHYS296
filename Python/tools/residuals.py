#!/usr/bin/env python
import sys,os
import numpy as np
import time
from numpy.random import choice,randn
from multiproc import MULTIPROC
import pandas as pd

class _RESIDUALS:

  # basic routimes to prepare the data sets

  def percent_to_absolute(self):
    for k in self.tabs:
      ucorr = [x for x in self.tabs[k] if '_u' in x and '%' in x]
      corr  = [x for x in self.tabs[k] if '_c' in x and '%' in x]
      if len(ucorr)!=0:
        for name in ucorr:
          mod_name=name.replace('%','')
          self.tabs[k][mod_name]=self.tabs[k]['value'] * self.tabs[k][name]
      if len(corr)!=0:
        for name in corr:
          mod_name=name.replace('%','')
          self.tabs[k][mod_name]=self.tabs[k]['value'] * self.tabs[k][name]/100.0

  def add_columns(self):
    for k in self.tabs:
      npts=len(self.tabs[k]['value'])
      self.tabs[k]['thy']=np.zeros(npts)
      self.tabs[k]['N']=np.zeros(npts)
      self.tabs[k]['residuals']=np.zeros(npts)
      self.tabs[k]['r-residuals']=np.zeros(npts)
      self.tabs[k]['Shift']=np.zeros(npts)

  def get_alpha(self):
    for k in self.tabs:
      npts=len(self.tabs[k]['value'])
      alpha2=np.zeros(npts) 
      ucorr = [x for x in self.tabs[k] if '_u' in x and '%' not in x]
      for kk in ucorr: 
        alpha2+=self.tabs[k][kk]**2
      
      self.tabs[k]['alpha']=alpha2**0.5

  def retrieve_norm_uncertainty(self):
    for k in self.tabs:
      norm  = [x for x in self.tabs[k] if '_c' in x and 'norm' in x]
      if len(norm)>1:
        raise ValueError('ERR: more than one normalization found at'%(k,kk))
      elif len(norm)==1:
        if self.conf['datasets'][self.reaction]['norm'][k]['fixed']==False:
          dN=self.tabs[k][norm[0]][0]/self.tabs[k]['value'][0]
          self.conf['datasets'][self.reaction]['norm'][k]['dN']=dN

  def setup_rparams(self):
    if 'rparams' not in self.conf: 
      self.conf['rparams']={}
    if self.reaction not in self.conf['rparams']: 
      self.conf['rparams'][self.reaction]={}
    for k in self.tabs:
      if k not in self.conf['rparams'][self.reaction]:
        self.conf['rparams'][self.reaction][k]={}
      corr = [x for x in self.tabs[k] if '_c' in x and '%' not in x]
      for c in corr: 
        self.conf['rparams'][self.reaction][k][c]={'value':0.0,'registered':False}

  def prepare_multiprocess(self):
    data=[]
    for k in self.tabs:
      for i in range(len(self.tabs[k]['value'])):
        data.append([k,i])
    if 'ncpus' in self.conf:
      ncpus=self.conf['ncpus']
    else: 
      ncpus=1
    self.mproc=MULTIPROC(ncpus,self._get_theory,data)

  # routines for IMC analysis
 
  def select_training_sets(self,tab):
    for k in tab:
      key=self.tabs[k].leys()[0]
      npts=len(self.tabs[k][key])
      tab[k]['iT']=np.zeros(npts)
      if npts>5:
        nptsT = int(self.conf['training frac']*npts)
        iT=choice(npts,nptsT,replace=False)
        for i in iT: tab[k]['iT'][i]=1       

  def resample(self,tab):
    for k in tab:
      key=self.tabs[k].leys()[0]
      npts=len(self.tabs[k][key])
      tab[k]+=randn(npts)*tab[k].alpha

  def setup_imc(self):
    # only useful for IMC
    if 'cross-validation' in self.conf:
      if self.conf['cross-validation']: self.select_training_sets()
    if 'bootstrap' in self.conf:
        if self.conf['bootstrap']: self.resample()

  # master setup

  def setup(self):
    self.percent_to_absolute()
    self.add_columns()
    self.get_alpha()
    self.retrieve_norm_uncertainty()
    self.setup_rparams()
    self.prepare_multiprocess()
    self.setup_imc()

  # residuals

  def _get_residuals(self,k):
    npts=len(self.tabs[k]['value'])
    exp=self.tabs[k]['value']
    norm=self.conf['datasets'][self.reaction]['norm'][k]['value']
    thy=self.tabs[k]['thy']/norm
    alpha=self.tabs[k]['alpha']
    corr = [x for x in self.tabs[k] if '_c' in x and '%' not in x and 'norm' not in x]
    N=np.ones(exp.size)
    ncorr=len(corr)
    if ncorr==0:
      self.tabs[k]['N']=N
      self.tabs[k]['residuals']=(exp-thy)/alpha
      self.tabs[k]['shift']=np.zeros(exp.size)
    else:
      B=[]
      beta=[]
      for c in corr:
        beta_=self.tabs[k][c] * (thy/exp)
        B.append(np.sum(beta_*(exp-thy)/alpha**2))
        beta.append(beta_)
      A=np.diag(np.diag(np.ones((ncorr,ncorr)))) \
        +np.einsum('ik,jk,k->ij',beta,beta,1/alpha**2)
      r=np.einsum('ij,j->i',np.linalg.inv(A),B)
      shift=np.einsum('k,ki->i',r,beta)
      for i in range(ncorr):
        self.conf['rparams'][self.reaction][k][corr[i]]['value']=r[i]
      self.tabs[k]['N']=N
      self.tabs[k]['residuals']=(exp-shift-thy)/alpha
      self.tabs[k]['shift']=shift
    return self.tabs[k]['residuals']

  def _get_rres(self,k):
    rres=[]
    rparams=self.conf['rparams'][self.reaction][k]
    for c in rparams:
      rres.append(rparams[c]['value'])
    return np.array(rres)

  def _get_nres(self,k):
    norm=self.conf['datasets'][self.reaction]['norm'][k]
    if 'dN' in norm:
      return (norm['value']-1)/norm['dN']
    else:
      return 0

  def get_theory(self):
    output=self.mproc.run()
    THY=[]
    for entry in output:
      k,i,thy=entry
      self.tabs[k]['thy'][i]=thy

  def get_npts(self):
    npts=0
    for k in self.tabs: 
      npts+=len(self.tabs[k]['value'])
    return npts

  def get_residuals(self):
    res,rres,nres=[],[],[]
    self.get_theory()
    for k in self.tabs: 
      res=np.append(res  ,self._get_residuals(k))
      rres=np.append(rres,self._get_rres(k))
      nres=np.append(nres,self._get_nres(k))
    return res,rres,nres

  # other functions

  def ___gen_report(self,verb=1,level=1):
    """
    verb = 0: Do not print on screen. Only return list of strings
    verv = 1: print on screen the report
    level= 0: only the total chi2s
    level= 1: include point by point 
    """

    L=[]

    L.append(self.reaction)

    for k in self.tabs:
      res =self._get_residuals(k)
      rres=self._get_rres(k)
      nres=self._get_nres(k)
      
      chi2=np.sum(res**2)
      rchi2=np.sum(rres**2)
      nchi2=nres**2
      tar=self.tabs[k]['target'].values[0]
      col=self.tabs[k]['col'].values[0].split()[0]
      npts=res.size
      L.append('%7d %10s %10s %5d %10.2f %10.2f %10.2f'%(k,tar,col,npts,chi2,rchi2,nchi2))

    if level==1:
      L.append('-'*100)  
      for k in self.conf['sidistab']:
        if len(self.conf['sidistab'][k].index)==0: continue 
        for i in range(len(self.conf['sidistab'][k].index)):
          x=self.conf['sidistab'][k]['x'].values[i]
          obs=self.conf['sidistab'][k]['obs'].values[i]
          Q2=self.conf['sidistab'][k]['Q2'].values[i]
          res=self.conf['sidistab'][k]['residuals'].values[i]
          thy=self.conf['sidistab'][k]['thy'].values[i]
          exp=self.conf['sidistab'][k]['value'].values[i]
          alpha=self.conf['sidistab'][k]['alpha'].values[i]
          rres=self.conf['sidistab'][k]['r-residuals'].values[i]
          col=self.conf['sidistab'][k]['col'].values[i]
          shift=self.conf['sidistab'][k].Shift.values[i]
          if res<0: chi2=-res**2
          else: chi2=res**2
          msg='%7s %7s x=%10.3e Q2=%10.3e exp=%10.3e alpha=%10.3e thy=%10.3e shift=%10.3e chi2=%10.3f'
          L.append(msg%(col,obs,x,Q2,exp,alpha,thy,shift,chi2))

    if verb==0:
      return L
    elif verb==1:
      for l in L: print l

  def ___save_results(self,path):
    save(self.tabs,path)


