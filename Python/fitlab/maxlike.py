#!/usr/bin/env python
import sys,os
import numpy as np
from tools.tools import load,save,checkdir,load_config
import time
from scipy.optimize import leastsq,minimize
import copy
from tools.config import conf
 
class ML:

  def __init__(self):
    self.cum_shifts=0
    self.iteration=0
    self.CHI2T=[]
    self.CHI2V=[]
    self.PARAMS=[]
    self.TOT=[]

    if 'screen mode' not in conf: 
      conf['screen mode']='plain'

  def get_stats(self,res,rres,nres,delay):

    shifts=conf['parman'].shifts
    etime = (time.time()-self.t0)/60

    npts=res.size
    chi2=np.sum(res**2)
    rchi2=np.sum(rres**2)
    nchi2=np.sum(nres**2)
    chi2tot=chi2+rchi2+nchi2
    dchi2=chi2tot-self.chi2tot
    if shifts>2: 
      if chi2tot<self.chi2tot:
        self.dchi2=self.chi2tot-chi2tot
        self.chi2tot=chi2tot


    self.status=[]
    self.status.append('JAM FITTER')
    self.status.append('count = %d'%self.cnt)
    self.status.append('elapsed time(mins)=%f'%etime)
    self.status.append('shifts  = %d'%shifts)
    self.status.append('npts    = %d'%npts)
    self.status.append('chi2    = %f'%chi2)
    self.status.append('rchi2   = %f'%rchi2)
    self.status.append('nchi2   = %f'%nchi2)
    self.status.append('chi2tot = %f'%(chi2tot))
    self.status.append('dchi2(iter)  = %f'%self.dchi2)
    self.status.append('dchi2(local) = %f'%dchi2)

    self.status.append('')
    self.status.extend(conf['resman'].gen_report())

    parstatus = conf['parman'].gen_report()

    if conf['screen mode']=='curses':

      self.status.append('')
      if    delay==False: self.status.append('')
      else: self.status.append('Its over!!!')
      self.status.append('')
      self.status.append('press q to exit')

      conf['screen'].clear()
      conf['screen'].border(0)
      if delay==False: conf['screen'].nodelay(1)
      else: conf['screen'].nodelay(0)

      for i in range(len(self.status)):
        conf['screen'].addstr(i+2,2,self.status[i])

      for i in range(len(parstatus)):
        conf['screen'].addstr(i+2,80,parstatus[i])


      conf['screen'].refresh()
      if conf['screen'].getch()==ord('q'):
        curses.endwin()
        self.gen_output()
        sys.exit()

    elif conf['screen mode']=='plain':
      for i in range(len(self.status)): print self.status[i]
      for i in range(len(parstatus)): print parstatus[i]
      if delay==True: 
        self.gen_output()
        sys.exit()

  def get_residuals(self,par,delay=False):
    self.cnt+=1
    res,rres,nres=conf['resman'].get_residuals(par)
    self.get_stats(res,rres,nres,delay)
    if len(rres)!=0: res=np.append(res,rres)
    if len(nres)!=0: res=np.append(res,nres)
    return res

  def gen_output(self):

    inputfile=conf['args'].config
    L=open(inputfile).readlines()
    for i in range(len(L)):

      if L[i].startswith('#'): continue

      if 'params' in L[i] and '<<' in L[i]:
        l=L[i].split('=')[0].replace('conf','').replace("'",'')
        k,kk=l.replace('][','@').replace('[','').replace(']','').split('@')[1:]
        left=L[i].split('<<')[0]
        right=L[i].split('>>')[1]
        L[i]=left+'<<%30.20e>>'%conf['params'][k.strip()][kk.strip()]['value']+right

      if 'norm' in L[i] and '<<' in L[i]:
        l=L[i].split('=')[0].replace('conf','').replace("'",'')
        dum1,k,dum2,kk=l.replace('][','@').replace('[','').replace(']','').split('@')
        left=L[i].split('<<')[0]
        right=L[i].split('>>')[1]

        value=conf['datasets'][k.strip()]['norm'][int(kk)]['value']
        L[i]=left+'<<%30.20e>>'%value+right

    #name=inputfile.split('/')[-1].replace('.py','')
    #outputdir='runs/%s'%name
    #checkdir(outputdir)
    #outputfile='%s/%s.py'%(outputdir,name)

    F=open(inputfile,'w')
    F.writelines(L)
    F.close()

  def get_chi2(self,par):
    return np.sum(self.get_residuals(par)**2)

  def run_minimize(self):

    guess=conf['parman'].par
    order=conf['parman'].order


    bounds=[]
    for entry in order:
      i,k,kk=entry
      if i==1:
        bounds.append([conf['params'][k][kk]['min'],conf['params'][k][kk]['max']])
      elif i==2:
        bounds.append([None,None])

    if conf['screen mode']=='curses':
      conf['screen']=curses.initscr()

    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0

    #self.get_residuals(guess)
    #sys.exit()

    res = minimize(self.get_chi2,guess, bounds=bounds,method='TNC')
    res=self.get_residuals(res.x,delay=True)

  def run_leastsq(self,gen_output=True):

    guess=conf['parman'].par
    order=conf['parman'].order

    bounds=[]
    for entry in order:
      i,k,kk=entry
      if i==1:
        bounds.append([conf['params'][k][kk]['min'],conf['params'][k][kk]['max']])
      elif i==2:
        bounds.append([None,None])

    if conf['screen mode']=='curses':
      conf['screen']=curses.initscr()

    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0
    fit=leastsq(self.get_residuals,guess,full_output = 1, ftol=1e-6)#,ftol=1e-2)#,factor=0.1)#,ftol=1e-2)
    res=self.get_residuals(fit[0])#,delay=True)
    if gen_output: self.gen_output()
    return fit[0]

  def analysis(self):
    self.gen_report()

  def gen_report(self):
    inputfile=conf['args'].config
    outdir='outputs/'+inputfile.split('/')[-1].replace('.py','')
    checkdir(outdir)
    par=conf['parman'].par
    self.t0 = time.time()
    self.cnt=0
    self.chi2tot=1e1000
    self.dchi2=0
    self.t0 = time.time()
    self.cnt=0
    self.get_residuals(par)
    report=conf['resman'].gen_report(verb=1,level=1)
    #save(report,'%s/report-ML'%outdir)

  def rap_fits(self):
    checkdir('maxlike')
    dy=[1.5,2.0,2.5,3.0,3.5]
    PAR={}
    for _dy in dy:
      print '#'*10
      conf['datasets']['sidis']['filters'][0]['filter']="z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9 and dy>%f"%_dy
      conf['resman'].setup()
      par=self.run_leastsq(gen_output=False)
      PAR[_dy]=par
    save(PAR,'maxlike/rap_fits.dat')

