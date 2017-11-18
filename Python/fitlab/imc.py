#!/usr/bin/env python
import sys,os
import numpy as np
from numpy.random import rand,randn,uniform
from scipy.integrate import quad,simps
from numpy import linalg as LA
from scipy.stats import gaussian_kde
import time
import cPickle 
import pickle 
import zlib
import warnings
warnings.filterwarnings('ignore')
from timeit import default_timer as timer
from numpy import linalg as la
from scipy.optimize import leastsq,minimize

# aux funcs

def lprint(msg):
  sys.stdout.write('\r')
  sys.stdout.write(msg)
  sys.stdout.flush()

def checkdir(path):
  if not os.path.exists(path): 
    os.makedirs(path)

def save(data,name):  
  compressed=zlib.compress(cPickle.dumps(data))
  f=open(name,"wb")
  f.writelines(compressed)
  f.close()

def load(name): 
  compressed=open(name,"rb").read()
  data=cPickle.loads(zlib.decompress(compressed))
  return data

def load_snapshot(name): 
  data=pickle.loads(open(name,"rb").read())
  self.__dict__.update(data.__dict__)

class ELLIPSE:

  def __init__(self,samples,kappa=1.0,N=None):

    self.N=N    
    self.dim=len(samples[0])
    # generate transformation matrix
    self.y0=np.mean(samples,axis=0)

    cov=self.get_cov(samples)
    w,v=np.linalg.eig(cov)
    icov=np.linalg.inv(cov)
    if np.any(np.isnan(icov)): raise ValueError('icov is nan')
    self.det=np.linalg.det(cov)

    v=np.transpose(v)
    for i in range(w.size): v[i]*=w[i]**0.5
    self.T=np.transpose(v)
    if np.any(np.isnan(self.T)): raise ValueError('T is nan')
    self.w=w
    self.v=v
    
    # get enlargement factor
    self.F=kappa*np.amax([np.einsum('i,ij,j',y-self.y0,icov,y-self.y0) for y in samples])**0.5
    if np.isnan(self.F): raise ValueError('F is nan')
    self.V=(self.F*self.det)**0.5
    self.gen_new_samples()

  def is_positive_semi_definite(self,cov):
    test=True
    w,v=np.linalg.eig(cov)
    if np.any(w<0): test=False
    if np.any(np.isnan(v)): test=False
    return test

  def get_cov(self,samples):
    cov=np.cov(np.transpose(samples))
    if self.is_positive_semi_definite(cov): 
      return cov
    else:
      sigma=np.abs(np.diagonal(cov))**0.5
      cnt=0
      while 1:
        cnt+=1
        if cnt%100==0: print 'fixing cov attempts:',cnt
        fake_samples=[sample+np.random.randn(sigma.size)*sigma for sample in samples]  
        cov=np.cov(np.transpose(fake_samples))
        if self.is_positive_semi_definite(cov): break
      return cov

  def gen_new_samples(self):

    # generate the unit sphere
    z=np.random.randn(self.N,self.dim)
    r=np.array([np.dot(z[i],z[i])**0.5 for i in range(self.N)])
    X=np.array([z[i]/r[i]*np.random.rand()**(1.0/self.dim) for i in range(self.N)])
    #print '='*100
    #print 'det='
    #print self.det
    #print 'w='
    #print self.w
    #print 'X='
    #print X[0]


    # generate sphere samples
    Y=np.einsum('ij,nj->ni',self.F*self.T,X) + self.y0
    #print Y[-10:]
    self.Y=[y for y in Y]

  def status(self):
    if len(self.Y)>0: return True
    else: return False

  def get_sample(self):
    return self.Y.pop()

class IMC:

  def __init__(self,conf):
    self.conf=conf
    self.get_nll = conf['nll']
    self.pmin=np.array([entry[0] for entry in conf['par lims']])
    self.pmax=np.array([entry[1] for entry in conf['par lims']])
    self.dp=np.array([entry[1]-entry[0] for entry in conf['par lims']])
    self.dim =len(conf['par lims']) 
    self.jac=np.prod(self.dp)
    self.N=conf['num points'] # number of active set
    self.factor=self.N/(self.N+1.) # factor to estimate prior mass
    self.V0=np.prod(self.dp)
    self.Vk=self.V0
    self.msg=''
    self.set_active_sets(self.N)
    self.status='ready'

  # param generators
  
  def gen_par(self):
    u=uniform(0,1,self.dim)
    return self.pmin + u*self.dp    

  def move(self,p0=None,verb=False):
    res = minimize(self.get_nll,p0, bounds=self.conf['par lims'],method='TNC')
    return res.x

  # sampling routine
  
  def set_active_sets(self,N):
    self.active_p=[]
    cnt_active=0
    while 1:
      lprint('getting initial active p: %d/%d'%(cnt_active+1,N))
      p=self.gen_par()
      cnt_active+=1
      self.active_p.append(p)
      if cnt_active==N: break

  def gen_sample(self):
    self.cnt+=1
    ellipse=ELLIPSE(self.active_p,self.conf['kappa'],self.N)
    self.Vk=ellipse.V
    for dum in range(len(self.active_p)):
      print dum
      self.active_p.pop(0)
      p=ellipse.get_sample().real
      pnew=self.move(p)
      self.active_p.append(pnew)

  def next(self,t_elapsed):
    self.gen_sample()
    msg='iter=%d  rel-err=%.3e  t-elapsed=%.3e    Vk/V0=%0.3e  '
    msg=msg%(self.cnt,t_elapsed,self.Vk/self.V0)
    lprint(msg)
    # stopping criterion
    if 'itmax' in self.conf and self.cnt==self.conf['itmax']: 
      self.status='stop'

  def results(self):
    data={}
    data['samples']=self.active_p
    data['weights']=np.ones(self.N)/float(self.N)
    return data

  def run(self):
    self.cnt=0  
    t1=timer()
    print 
    while 1:
      try:
        t2=timer()
        self.next(t2-t1)
        if self.status=='stop': break
      except KeyboardInterrupt:
        break
    print 
    return self.results() 

def example1():

  dim=2
  sigma=np.ones(dim)
  mean=np.zeros(dim)*0.1
  def likelihood(p):
    norm=1/np.prod(np.sqrt(2*np.pi*sigma**2))
    return norm * np.exp(-0.5*np.sum((p-mean)**2/sigma**2))

  nll=lambda p: -np.log(likelihood(p))


  print 'True:'
  print 'logz=',np.log(1)
  print 'min nll=',nll(np.zeros(dim))
  print 
  #print 'VEGAS:'
  #import vegas
  #integ = vegas.Integrator([[-5, 5] for i in range(dim)])
  #result = integ(likelihood, nitn=10, neval=1000)
  #print result.summary() 
  #print 'logz = %s    Q = %.2f' % (np.log(result), result.Q)

  print 
  print 'Nested Sampling:'
  conf={}
  conf['nll'] = nll
  conf['par lims'] =[[-5,5] for i in range(dim)]
  conf['tol']=1e-10
  conf['num points'] = 50

  conf['method']='cov'
  conf['kappa']=1.0
  conf['sample size']= 100

  #conf['method']='kde'
  #conf['kde bw']=None

  NEST(conf).run()

if __name__=='__main__':

  example1()







