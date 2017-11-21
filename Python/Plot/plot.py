#!/usr/bin/env python
"""
Created on Tue Nov 21 09:23:42 2017

@author: JohnJohn
"""
import os
#from openpyxl import Workbook as wb
#from xlrd import open_workbook as wb 
from pandas import read_excel as re
import matplotlib.pyplot as plt

def delta(stat_u,sys_u):
    error = (stat_u**2.0 + sys_u**2.0)**(1/2.0)
    return error



#Path = 'C:\Users\JohnJohn\Documents\Backup\School\Fa_17\Phys_296\PSU-PHYS296\Python\database\sidis\expdata'
#F=os.listdir('C:\Users\JohnJohn\Documents\Backup\School\Fa_17\Phys_296\PSU-PHYS296\Python\database\sidis\expdata')
#f = Path+F[-1]

Path = 'C:/Users/JohnJohn/Documents/Backup/School/Fa_17/Phys_296/PSU-PHYS296/Python/database/sidis/expdata/2001.xlsx'
#workbook = wb(Path)
#sheet = workbook.sheet_by_index(0)
data = re(Path)
delta = delta(data.stat_u,data.sys_u)
plt.errorbar(data.pT,data.value,delta,fmt='o',mfc='green',ecolor='green',capsize=5) # maker vs fmt
plt.title('Figure 1')
plt.ylabel('ALL')
plt.xlabel('pT')

