#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 15:43:22 2017

@author: avp5627
"""

#!/usr/bin/env python
import os
import pandas as pd

F=os.listdir('./')
F=[f for f in F if f.endswith('.dat')]

TABLE=[]
for f in F:
  L=open(f).readlines()
  L=[l.strip() for l in L]
  H=L[1].replace('#PhT','PhT').split()
  T=L[2:]
  T=[[float(x) for x in l.split()] for l in T if l!='']
  TABLE.extend(T)

PD=pd.DataFrame(TABLE,columns=H)
print PD

writer = pd.ExcelWriter('expdata.xlsx')
PD.to_excel(writer, index=False,header=True)
