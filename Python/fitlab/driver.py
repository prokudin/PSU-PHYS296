#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import json


data={}
data['model']="WW+Gaussian"
data['particle']="pi+"
data['target']="proton"
data['varialves']=["Fuu"]
data['axis']=[]
data['axis'].append({ "name": "a", "bins":  40, "min":  0.025, "max":   0.995, "scale":"arb" ,"description":"Bjorken x"})
data['axis'].append({ "name": "b", "bins":  40, "min":  0.95,  "max":    20.0, "scale":"arb", "description":"Q^2"},)
data['axis'].append({ "name": "c", "bins":  40, "min":  0.025, "max":   0.995, "scale":"lin", "description":"hadron frac. energy z"},)
data['axis'].append({ "name": "d", "bins":  40, "min":  0.00,  "max":    2.00, "scale":"lin", "description":"transverse momentum PT"})
data = json.dumps(data)
f = open("test.json",'w')
f.writelines(data)
f.close()

