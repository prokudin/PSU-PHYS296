#!/usr/bin/env python
import sys,os
import numpy as np
import time

def get_cov(sample):
  n=len(sample)
  mean=np.mean(sample,axis=0)
  diff=[s-mean for s in sample]
  return np.einsum('ki,kj->ij',diff,diff)/(n-1)


@profile
def test(samples):                                
  cov2=get_cov(samples)
  cov1=np.cov(np.transpose(samples))
  #del cov1
  time.sleep(10)
  return None

samples=np.random.rand(50,20)
test(samples)                           




