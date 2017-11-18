#!/usr/bin/env python
import sys,os
import numpy as np
from timeit import default_timer as timer
from tools.tools import lprint

class SPEEDTEST:

  def __init__(self,conf):
    self.conf=conf

  def run(self):
    resman=self.conf['resman']
    parman=self.conf['parman']
    args=self.conf['args']

    a = timer()
    np.random.seed(12345)
    for i in range(10):
      msg='[%d/10]'%i
      lprint(msg)
      R=np.random.randn(len(parman.par))
      par=parman.par+1e-3*R
      res,rres,nres=resman.get_residuals(par)
    print '\nt-elapsed (sec): ',(timer()-a)/10
    print 'npts=%d'%res.size
    print 'chi2=%f'%np.sum(res**2)
    print 'chi2/npts=%f'%(np.sum(res**2)/res.size)
    print 'chi2(r)=%f'%np.sum(rres**2)
    print 'chi2(norm)=%f'%np.sum(nres**2)

