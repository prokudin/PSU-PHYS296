#!/usr/bin/env python
import sys,os
import numpy as np
import pandas as pd
from tools import isnumeric

class _READER:

  def __init__(self,conf):
    self.conf=conf

  def apply_cuts(self,tab,k):
    if 'filters' in self.conf['datasets'][self.reaction]:
      filters=self.conf['datasets'][self.reaction]['filters']
      for i in filters:
        if k in filters[i]['idx']:
          tab=tab.query(filters[i]['filter'])
    return tab

  def load_data_sets(self,reaction):
    self.reaction=reaction
    if reaction not in self.conf['datasets']: return None
    XLSX=self.conf['datasets'][reaction]['xlsx']
    TAB={}
    for k in XLSX: 
      sys.stdout.write('\r')
      sys.stdout.write('loading %s data sets %d'%(reaction,k))
      sys.stdout.flush()
      fname=self.conf['datasets'][reaction]['xlsx'][k]
      tab=pd.read_excel(fname)
      tab=self.modify_table(tab,k)
      TAB[k]=tab.to_dict(orient='list')
      for kk in TAB[k]:
        if len(TAB[k][kk])==0: continue
        if isnumeric(TAB[k][kk][0]):
          TAB[k][kk]=np.array(TAB[k][kk])

    return TAB

