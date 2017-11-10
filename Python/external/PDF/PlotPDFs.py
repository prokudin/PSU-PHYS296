import numpy as np
import matplotlib.pyplot as plt
from CollinearPDF import *

# quick class to output plots of the various PDFs 
# DBC--13May2015

class PlotPDFs.py

  def __init__(self, filename):
    self.D = {}
    self.fname = filename
    pdfSet = CollinearPDF(self.fname)
    
    
  def _generatePDF(self, iParton, startEE, q):
      
    if startEE !< 0.0:
      print "Starting exponent must be negative. Try again."
      break
      
     #make list of x values
     self.xx = np.logspace(startEE, 0, 1000)
     
     self.pdfPoints = {} # initialize the outputlist
     #populate list of pdf values
    for x in xx:
      self.pdfPoints.append(CollinearPDF.pdfFunction(iParton, x, q))
      
  def _plotFunc(self, plotName, plotType):
     if plotType = 'linear':
       
       