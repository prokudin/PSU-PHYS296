#!/usr/bin/env python
import sys,os
import numpy as np
#from mpmath import fp,mp
import time

def get_cov(sample):
  return np.cov(sample.T)
  #n=len(sample)
  #mean=np.mean(sample,axis=0)
  #diff=[s-mean for s in sample]
  #return np.einsum('ki,kj->ij',diff,diff)/(n-1)


@profile
def test():
  for i in range(10):
    samples=np.random.rand(50,3)
    cov=np.cov(samples.T)
  #w,v=np.linalg.eig(cov)
  #time.sleep(10)
  #w, v = mp.eig(mp.matrix(cov))
  #w=[float(x) for x in w]
  #v=np.array(v.tolist(),dtype=float)
  #print np.einsum('li,ij,jm->lm',v.T,cov,v)
  #print w

test()
#test()


