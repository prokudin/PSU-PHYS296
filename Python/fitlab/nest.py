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

  def __init__(self,samples,kappa=1.0,iteration=None,N=None):
    self.iteration=iteration
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

  # fix cov matrix. 

  def fix_cov1(self,samples,cov):
    sigma=np.abs(np.diagonal(cov))**0.5
    cnt=0
    while 1:
      cnt+=1
      #if cnt%100==0: 
      print '\nfixing cov attempts:',cnt
      fake_samples=[sample+np.random.randn(sigma.size)*sigma for sample in samples]  
      cov=np.cov(np.transpose(fake_samples))
      if self.is_positive_semi_definite(cov): break
    return cov

  def vol_prefactor(self,n):
      """Volume constant for an n-dimensional sphere:
      for n even:      (2pi)^(n    /2) / (2 * 4 * ... * n)
      for n odd :  2 * (2pi)^((n-1)/2) / (1 * 3 * ... * n)
      """
      if n % 2 == 0:
          f = 1.
          i = 2
          while i <= n:
              f *= (2. / i * np.pi)
              i += 2
      else:
          f = 2.
          i = 3
          while i <= n:
              f *= (2. / i * np.pi)
              i += 2
      return f

  def make_eigvals_positive(self,cov, targetprod):
      """
      For the symmetric square matrix ``cov``, increase any zero eigenvalues
      to fulfill the given target product of eigenvalues.
      Returns a (possibly) new matrix.
      """  
      w, v = np.linalg.eigh(cov)  # Use eigh because we assume a is symmetric.
      mask = w < 1.e-10
      if np.any(mask):
          nzprod = np.product(w[~mask])  # product of nonzero eigenvalues
          nzeros = mask.sum()  # number of zero eigenvalues
          w[mask] = (targetprod / nzprod) ** (1./nzeros)  # adjust zero eigvals
          cov = np.dot(np.dot(v, np.diag(w)), np.linalg.inv(v))  # re-form cov
      return cov

  def fix_cov2(self,samples,cov):
      """
      Ensure that ``cov`` is nonsingular.
      It can be singular when the ellipsoid has zero volume, which happens
      when npoints <= ndim or when enough points are linear combinations
      of other points. (e.g., npoints = ndim+1 but one point is a linear
      combination of others). When this happens, we expand the ellipse
      in the zero dimensions to fulfill the volume expected from
      ``pointvol``.
      """
      print '\nfixing cov'
      npoints=self.N
      expected_vol = np.exp(-self.iteration / float(npoints) )
      pointvol= expected_vol / npoints 
      targetprod = (npoints * pointvol / self.vol_prefactor(self.dim))**2
      return self.make_eigvals_positive(cov, targetprod)

  def get_cov(self,samples):
    cov=np.cov(np.transpose(samples))
    if self.is_positive_semi_definite(cov): 
      return cov
    else:
      return self.fix_cov1(samples,cov)
      #return self.fix_cov2(samples,cov)

      print 
      print 'cov is singular'
      sys.exit() 

  def gen_new_samples(self):

    # generate the unit sphere
    z=np.random.randn(self.N,self.dim)
    r=np.array([np.dot(z[i],z[i])**0.5 for i in range(self.N)])
    X=np.array([z[i]/r[i]*np.random.rand()**(1.0/self.dim) for i in range(self.N)])

    # generate sphere samples
    Y=np.einsum('ij,nj->ni',self.F*self.T,X) + self.y0
    #print Y[-10:]
    self.Y=[y for y in Y]

  def status(self):
    if len(self.Y)>0: return True
    else: return False

  def get_sample(self):
    return self.Y.pop()

class NEST:

  def __init__(self,conf):
    self.conf=conf
    self.get_nll = conf['nll']
    self.pmin=np.array([entry[0] for entry in conf['par lims']])
    self.pmax=np.array([entry[1] for entry in conf['par lims']])
    self.dp=np.array([float(entry[1]-entry[0]) for entry in conf['par lims']])
    self.dim =len(conf['par lims']) 
    self.jac=np.abs(np.prod(self.dp))
    self.N=conf['num points'] # number of active set
    self.factor=self.N/(self.N+1.) # factor to estimate prior mass
    self.V0=np.prod(self.dp)
    self.Vk=self.V0
    self.msg=''

    if 'nestout' not in conf:
      self.samples_p=[]
      self.samples_x=[1]
      self.samples_l=[0]
      self.samples_nll=[]
      self.logz=[]
      self.cnt=0 # counter
      self.attempts1=None
      self.attempts2=None
      self.set_active_sets(self.N)
    else:
      self.active_p=conf['nestout']['active p']
      self.active_nll=conf['nestout']['active nll']
      self.samples_p=conf['nestout']['samples'][::-1]
      self.samples_x=conf['nestout']['x'][::-1]
      self.samples_l=conf['nestout']['l'][::-1]
      self.samples_nll=[-np.log(l) for l in self.samples_l]
      self.logz=conf['nestout']['logz']
      self.cnt=len(self.samples_l) # counter
      self.attempts1=None
      self.attempts2=None

    self.status='ready'

  # param generators
  
  def gen_par(self):
    u=uniform(0,1,self.dim)
    return self.pmin + u*self.dp    
  
  def gen_par_flat(self,nll,verb=False):
    self.attempts1=0
    self.attempts2=0
    pmin=np.amin(self.active_p,axis=0)
    pmax=np.amax(self.active_p,axis=0)
    dp=(pmax-pmin)
    self.Vk=np.prod(dp)
    dp=(pmax-pmin)*self.conf['kappa']
    pmid=0.5*(pmin+pmax)
    pmin=pmid-dp/2
    while 1:
      self.attempts2+=1

      if self.attempts2>1000:
        pmin=np.amin(self.active_p,axis=0)
        pmax=np.amax(self.active_p,axis=0)
        dp=(pmax-pmin)
        self.conf['kappa']=1

      u=uniform(0,1,self.dim)
      p=pmin + u*dp
      if any([x<0 for x in p-self.pmin]) or any([x<0 for x in self.pmax-p]): continue
      _nll = self.get_nll(p)
      if _nll<nll: break
    return p,_nll 
  
  def gen_par_kde(self,nll):
    kde=gaussian_kde(np.transpose(self.active_p),self.conf['kde bw'])
    self.attempts=0
    while 1:
      p=np.transpose(kde.resample(1))[0]
      if any([x<0 for x in p-self.pmin]): continue
      if any([x<0 for x in self.pmax-p]): continue    
      self.attempts+=1
      _nll = self.get_nll(p)
      if _nll<nll: break
    return p,_nll
  
  def gen_par_hmc(self,par,nll):
    dim=len(par)  
    delta=0.1#self.conf['hmc delta']
    steps=30#self.conf['hmc steps']
  
    get_U=lambda q: self.get_nll(q)
    
    def get_dU_i(i,q):
      h=0.01*q
      shift=np.zeros(dim)
      shift[i]=h[i]
      return (get_U(q+shift)-get_U(q-shift))/2/h[i]
    
    get_dU=lambda q: np.array([get_dU_i(i,q) for i in range(dim)])
    get_K=lambda p: np.dot(p,p)/2
    get_dK=lambda p: p
    get_H=lambda q,p: get_U(q) + get_K(p)
  
    q0=np.copy(par)
    H0=get_H(q0,par)
    while 1:
      print 'hmc attempt',self.attempts
      p0=randn(dim)
      p=np.copy(p0)
      q=np.copy(q0)
      p=p0-delta/2*get_dU(q0)
      q=q0+delta*p
      for i in range(steps):
        print 'walking',i
        p=p-delta*get_dU(q)
        q=q+delta*get_dK(p)
      p=p-delta/2*get_dU(q)
      _nll = self.get_nll(q)
      if _nll<nll: break
    return q,_nll

  def gen_par_cov(self,nll,p0=None,verb=False):
    self.attempts1=0
    self.attempts2=0
    ellipse=ELLIPSE(self.active_p,self.conf['kappa'],self.cnt,self.conf['sample size'])
    self.Vk=ellipse.V
    cnt=0
    while 1:
      if ellipse.status()==False: 
        ellipse.gen_new_samples()
      p=ellipse.get_sample().real
      if np.any(np.isnan(p)):
        raise ValueError('parameters are nan')
      self.attempts1+=1
      if any([x<0 for x in p-self.pmin]) or any([x<0 for x in self.pmax-p]): continue
      self.attempts2+=1
      _nll = self.get_nll(p)
      if _nll<nll: break
    return p,_nll

  # nested sampling routines
  
  def gen_sample(self):
  
    # remove entry from active arrays
    imax=np.argmax(self.active_nll)
    nll=self.active_nll.pop(imax)
    p=self.active_p.pop(imax)

    # update samples 
    self.samples_nll.append(nll)
    self.samples_l.append(np.exp(-nll))
    if len(self.samples_p)==0:
      self.samples_p=[p]
    else:
      self.samples_p=np.concatenate((self.samples_p,[p]),axis=0)
    self.samples_x.append(self.samples_x[-1]*self.factor)  
    self.cnt+=1
  
    _p,_nll=self.gen_par_flat(nll)
    #if self.cnt<self.conf['burn size']:
    #  _p,_nll=self.gen_par_flat(nll)
    #else:
    #  _p,_nll=self.gen_par_cov(nll,p)
  
    self.active_nll.append(_nll)
    self.active_p.append(_p)
  
  def get_logz(self):
    return np.log(np.trapz(self.samples_l[::-1],self.samples_x[::-1])*self.jac)
  
  def set_active_sets(self,N):
    self.active_p=[]
    self.active_nll=[]
    cnt_active=0
    while 1:
      lprint('getting initial active p: %d/%d'%(cnt_active+1,N))
      p=self.gen_par()
      nll=self.get_nll(p)
      cnt_active+=1
      self.active_p.append(p)
      self.active_nll.append(nll)
      if cnt_active==N: break

  def next(self,t_elapsed):
    self.gen_sample()
    self.logz.append(self.get_logz())
    if self.cnt>2: 
      #z_past=np.average([np.exp(self.logz[-2]),np.exp(self.logz[-3]),np.exp(self.logz[-4])])
      std=np.std(np.array(self.active_nll)+self.conf['data size'])
      mean=np.mean(np.array(self.active_nll)+self.conf['data size'])
      #z_current=np.exp(self.logz[-1])
      #rel = np.abs(1-z_past/z_current)
      rel = std/mean
      nllmax=np.amax(self.active_nll)
      nllmin=np.amin(self.active_nll)
      msg='iter=%d  logz=%.3f rel-err=%.3e  t-elapsed=%.3e  nll_min=%.3e nll_max=%0.3e  attemps1=%10d  attemps2=%10d Vk/V0=%0.3e  %s  '
      msg=msg%(self.cnt,self.logz[-1],rel,t_elapsed,nllmin,nllmax,self.attempts1,self.attempts2,self.Vk/self.V0,self.msg)
      lprint(msg)
      # stopping criterion
      if 'itmax' in self.conf and self.cnt==self.conf['itmax']: 
        self.status='stop'
      if self.logz[-1]>-1e10  and self.conf['tol']!=None:
        if rel<self.conf['tol']: self.status='stop'
        #l=np.exp(-self.active_nll[0])
        #x=self.samples_x[-1]
        #dz=l*x
        #if np.log(dz)<self.logz[-1]+np.log(self.conf['tol']):
        #  self.status='stop'

  def results(self):
    dx=0.5*(np.array(self.samples_x[:-1])-np.array(self.samples_x[1:]))
    log_l=np.log(self.samples_l[1:])
    log_w=np.log(dx*self.jac)+log_l
    weights=np.exp(log_w)
    weights/=np.sum(weights)
  
    data={}
    data['samples']=self.samples_p[::-1]
    data['x']=self.samples_x[1:][::-1]
    data['l']=self.samples_l[1:][::-1]
    data['logz']=self.logz
    data['weights']=weights[::-1]
    data['active p']=self.active_p
    data['active nll']=self.active_nll
    return data

  def run(self):
  
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

  dim=5
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

def example2():

  shift1=np.array([1, 1])
  shift2=np.array([1,-1])
  widths=[1,1,1,1]
  def likelihood(a):
    return 1./(2*np.pi*widths[0]**2)* np.exp(-np.sum((a-shift1)**2/2/widths[0]**2))\
         + 1./(2*np.pi*widths[1]**2)* np.exp(-np.sum((a+shift1)**2/2/widths[1]**2))\
         + 1./(2*np.pi*widths[2]**2)* np.exp(-np.sum((a-shift2)**2/2/widths[2]**2))\
         + 1./(2*np.pi*widths[3]**2)* np.exp(-np.sum((a+shift2)**2/2/widths[3]**2))
  nll=lambda p: -np.log(likelihood(p))


  print 'True:'
  print 'logz=',np.log(4)


  print 
  print 'Nested Sampling:'
  conf={}
  conf['nll'] = nll
  conf['par lims'] =[[-20,20],[-20,20]]
  conf['tol']=1e-4
  conf['num points'] = 100


  conf['method']='cov'
  conf['kappa']=1

  #conf['method']='kde'
  #conf['kde bw']=None
  ##conf['itmax']=None

  conf['sample size']= N
  print 'N=',N
  NEST(conf).run()

  #print 
  #print 'VEGAS:'
  #import vegas
  #integ = vegas.Integrator([[-20, 20], [-20, 20]])
  #result = integ(likelihood, nitn=10, neval=1000)
  #print result.summary() 
  #print 'logz = %s    Q = %.2f' % (np.log(result), result.Q)

if __name__=='__main__':

  example1()







