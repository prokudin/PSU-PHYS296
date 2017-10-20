#!/usr/bin/env python
import sys,os
import numpy as np
import pandas as pd
from tools.reader import _READER

class READER(_READER):

  def __init__(self,conf):
    self.conf=conf

  def modify_table(self,tab,k):
    tab=self.apply_cuts(tab,k)
    return tab

if __name__ == "__main__":

  conf={}
  conf['datasets']={}
  conf['datasets']['sidis']={}

  conf['datasets']['sidis']['xlsx']={}
  conf['datasets']['sidis']['xlsx'][5000]='../../database/sidis/expdata/5000.xlsx' 
  conf['datasets']['sidis']['xlsx'][5008]='../../database/sidis/expdata/5008.xlsx' 


  conf['datasets']['sidis']['filters']={}
  conf['datasets']['sidis']['filters'][1]={}
  conf['datasets']['sidis']['filters'][1]['list']=range(1000,2000)
  conf['datasets']['sidis']['filters'][1]['cond']=[]
  conf['datasets']['sidis']['filters'][1]['cond'].append("z<0.6")
  conf['datasets']['sidis']['filters'][1]['cond'].append("Q2>1.69")
  conf['datasets']['sidis']['filters'][1]['cond'].append("pT>0.2 and pT<0.9")

  TAB=READER(conf).load_data_sets('sidis')
  print TAB





