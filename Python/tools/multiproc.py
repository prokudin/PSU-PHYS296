#!/usr/bin/env python
import sys,os
from multiprocessing import Process,Queue,Pool,Pipe

class MULTIPROC:

  def __init__(self,ncores,func,data):
    self.ncores=ncores
    self.func=func
    self.data=data
    n=len(self.data)/(self.ncores)
    self.ini,self.fin=[],[]
    for i in range(self.ncores):
      self.ini.append(i*n)
      self.fin.append((i+1)*n)
    self.fin[-1]=len(self.data)
    self.pipes=[Pipe() for i in range(self.ncores)]

  def func_wrap(self,data,pipe):
    output=[self.func(entry) for entry in data]
    pipe.send(output)

  def singlecore(self):
    return [self.func(entry) for entry in self.data]

  def multicore(self):
    P = [Process(target=self.func_wrap, \
                      args=(self.data[self.ini[i]:self.fin[i]],self.pipes[i][1],)) \
                      for i in range(self.ncores)]
    for p in P: p.start() 
    for p in P: p.join() 
    output=[]
    for pipe in self.pipes: output.extend(pipe[0].recv())
    return output

  def run(self):
    if self.ncores==1: 
      return self.singlecore()
    else:
      return self.multicore()

