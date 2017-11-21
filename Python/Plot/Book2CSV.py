#!/usr/bin/env python
"""
Created on Tue Nov 21 11:23:11 2017

@author: JohnJohn
"""
import os
from pandas import read_excel, DataFrame as re, df


def Book2CSV(WorkBook_Path):
    dat = re(WorkBook_Path) # Read .xlsx
    data = df(dat)  # panda.dataFrame
    data.to_csv("class_data_2.csv",index=False)    # save as csv with out index vaues
    dir = os.getcwd()
    return print()
    
