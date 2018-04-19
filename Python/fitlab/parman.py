#!/usr/bin/env python
import sys,os
import numpy as np 
from numpy.random import choice,randn,uniform
#from tools.tools import load_config
from tools.config import conf,load_config
import pandas as pd
from external.PDF.CT10 import CT10
#from external.CJLIB.CJ import CJ
from external.DSSLIB.DSS import DSS
from external.LSSLIB.LSS import LSS
from qcdlib.tmdlib import PDF,PPDF,FF
from qcdlib.aux import AUX
from qcdlib.alphaS import ALPHAS
from obslib.dis.stfuncs import STFUNCS as DIS_STFUNCS
from obslib.sidis.reader import READER as SIDIS_READER
from obslib.sidis.stfuncs import STFUNCS as SIDIS_STFUNCS
from obslib.sidis.residuals import RESIDUALS as SIDIS_RESIDUALS

class PARMAN:

  def __init__(self):
    self.get_ordered_free_params()
    self.set_new_params(self.par,initial=True)

  def get_ordered_free_params(self):
    self.par=[]
    self.order=[]

    for k in conf['params']:
      for kk in conf['params'][k]:
        if conf['params'][k][kk]['fixed']==False:
          #p=uniform(conf['params'][k][kk]['min'],conf['params'][k][kk]['max'],1)[0]
          self.par.append(conf['params'][k][kk]['value'])
          self.order.append([1,k,kk])

    for k in conf['datasets']:
      for kk in conf['datasets'][k]['norm']:
        if conf['datasets'][k]['norm'][kk]['fixed']==False:
          self.par.append(conf['datasets'][k]['norm'][kk]['value'])
          self.order.append([2,k,kk])
        
  def set_new_params(self,parnew,initial=False):
    self.shifts=0
    semaphore={}

    for i in range(len(self.order)):
      ii,k,kk=self.order[i]  
      if ii==1:
        if k not in semaphore: semaphore[k]=0
        if conf['params'][k][kk]['value']!=parnew[i]:
          conf['params'][k][kk]['value']=parnew[i]
          semaphore[k]=1
          self.shifts+=1
      elif ii==2:
        if conf['datasets'][k]['norm'][kk]['value']!=parnew[i]:
          conf['datasets'][k]['norm'][kk]['value']=parnew[i]
          self.shifts+=1

    if initial:
      #for k in semaphore: semaphore[k]=1
      for k in conf['params']: semaphore[k]=1

    self.propagate_params(semaphore)

  def gen_report(self):
    L=[]

    for k in conf['params']:
      for kk in sorted(conf['params'][k]):
        if conf['params'][k][kk]['fixed']==False: 
          if conf['params'][k][kk]['value']<0:
            L.append('%-10s  %-20s  %10.5e'%(k,kk,conf['params'][k][kk]['value']))
          else:
            L.append('%-10s  %-20s   %10.5e'%(k,kk,conf['params'][k][kk]['value']))

    for k in conf['datasets']:
      for kk in conf['datasets'][k]['norm']:
        if conf['datasets'][k]['norm'][kk]['fixed']==False: 
          L.append('%10s %10s %10d  %10.5e'%('norm',k,kk,conf['datasets'][k]['norm'][kk]['value']))
    return L

  def propagate_params(self,semaphore):
    #print 'semaphore:',semaphore
    if 'pdf' in semaphore and semaphore['pdf']==1: self.set_pdf_params()
    if 'ff'  in semaphore and semaphore['ff']==1:  self.set_ff_params()
    if 'ppdf'  in semaphore and semaphore['ppdf']==1:  self.set_ppdf_params()
    if 'gk'  in semaphore and semaphore['gk']==1:  self.set_gk_params()

  def set_ppdf_params(self):
    conf['ppdf'].widths0['valence']=conf['params']['ppdf']['widths0 valence']['value']
    conf['ppdf'].widths0['sea']=conf['params']['ppdf']['widths0 sea']['value']
    conf['ppdf'].setup() 

  def set_pdf_params(self):
    conf['pdf'].widths0['valence']=conf['params']['pdf']['widths0 valence']['value']
    conf['pdf'].widths0['sea']=conf['params']['pdf']['widths0 sea']['value']
    conf['pdf'].setup() 
  
  def set_ff_params(self):
    conf['ff'].widths0['pi+ fav']=conf['params']['ff']['widths0 pi+ fav']['value']
    conf['ff'].widths0['pi+ unfav']=conf['params']['ff']['widths0 pi+ unfav']['value']
    conf['ff'].widths0['k+ fav']=conf['params']['ff']['widths0 k+ fav']['value']
    conf['ff'].widths0['k+ unfav']=conf['params']['ff']['widths0 k+ unfav']['value']
    conf['ff'].setup() 

  def set_gk_params(self):
    conf['gk'].Q0=conf['params']['gk']['Q0']['value']
    conf['gk'].g2=conf['params']['gk']['g2']['value']
    conf['gk'].bmax=conf['params']['gk']['bmax']['value']
    conf['gk'].bmin=conf['params']['gk']['bmin']['value']
  


if __name__=='__main__':

  conf=load_config('input.py')

  # setup for inclusive dis
  conf['alphaSmode']='backward'
  conf['Q20']=1
  conf['order']='NLO'
  conf['aux']=AUX()
  conf['alphaS']=ALPHAS()
  conf['pdf-NLO']=CT10()
  conf['dis stfuncs']=DIS_STFUNCS()

  # setup tmds
  conf['order']='LO'
  conf['_pdf']=CJ()
  conf['_ppdf']=LSS()
  conf['_ff']=DSS()
  conf['pdf']=PDF()
  conf['ppdf']=PPDF()
  conf['ff']=FF()
 
  # setup sidis
  conf['sidis tabs']=SIDIS_READER().load_data_sets('sidis')
  conf['sidis stfuncs']=SIDIS_STFUNCS()
  conf['sidis residuals']=SIDIS_RESIDUALS()

  parman=PARMAN()




