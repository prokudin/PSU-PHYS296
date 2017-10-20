#!/usr/bin/env python
import sys,os
from tools import load

class STORAGE(object):

  def __init__(self,path,fname):
    checkdir(path)
    self.fname=path+'/'+fname
    self.load_storage(self.fname)

  def load_storage(self,fname):

    # create storage if not created
    if not os.path.exists(fname):
      #print fname+' is been created !'
      self.storage={}
    else:
      self.storage=load(fname)

  def gen_key(self,x,Q2):
    return 'x='+str(x)+',Q2='+str(Q2)

  def query(self,x,Q2):
    key=self.gen_key(x,Q2)
    return any([key==k for k in self.storage.keys()])

  def register(self,func,x,Q2,name):

    if self.query(x,Q2)==False:
      #print 'register x=%0.3e Q2=%0.3e at %s'%(x,Q2,name)
      key=self.gen_key(x,Q2)
      self.storage[key]=func(x,Q2)

  def Save(self):
    save(self.storage,self.fname)

  def retrieve(self,x,Q2):
    key=self.gen_key(x,Q2)
    return self.storage[key]


