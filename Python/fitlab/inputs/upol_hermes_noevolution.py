conf={}

############################################################################
# resouce allocation

conf['ncpus']=1

############################################################################
# maxlike setup

conf['screen mode']='plain'
#conf['screen mode']='curses'

############################################################################
# mcsetup

conf['method']='kde'
conf['itmax']=None
conf['tol']=1e-6
conf['kde bw']=None
conf['num points factor']=10


############################################################################
# paths to external

#conf['path2CJ'] ='../external/CJLIB'
conf['path2CT10'] ='../external/PDF'
conf['path2LSS']='../external/LSSLIB'
conf['path2DSS']='../external/DSSLIB'

############################################################################
# params

conf['params']={}

conf['params']['pdf']={}
conf['params']['pdf']['widths0 valence']  = {'value':<<    2.62111526535736771848e-01>>,'fixed':False,'min':0,'max':10}
conf['params']['pdf']['widths0 sea']      = {'value':<<    4.54546309613160559593e-01>>,'fixed':False,'min':0,'max':10}

conf['params']['ff']={}
conf['params']['ff']['widths0 pi+ fav']   = {'value':<<    2.11195846199639658547e-01>>,'fixed':False,'min':0,'max':1}
conf['params']['ff']['widths0 pi+ unfav'] = {'value':<<    1.93054237855291910275e-01>>,'fixed':False,'min':0,'max':1}
conf['params']['ff']['widths0 k+ fav']    = {'value':<<    2.98635175670050156960e-01>>,'fixed':False,'min':0,'max':1}
conf['params']['ff']['widths0 k+ unfav']  = {'value':<<    1.38449119642011064801e-01>>,'fixed':False,'min':0,'max':1}

conf['params']['gk']={}
conf['params']['gk']['g2']   = {'value':<<    0.00000000000000000000e+00>>,'fixed':False,'min':0,'max':10}
conf['params']['gk']['Q0'] = {'value':<<    1.00000000000000000000e+00>>,'fixed':True,'min':1,'max':10}
conf['params']['gk']['bmax'] = {'value':<<    1.30756493325539491224e+00>>,'fixed':False,'min':0.1,'max':10}#0.01,'max':10}
conf['params']['gk']['bmin'] = {'value':<<    3.95183715078905706264e-01>>,'fixed':True,'min':0.0001,'max':10}#0.01,'max':10}


############################################################################
# set data sets

conf['datasets']={}

# SIDIS

conf['datasets']['sidis']={}


conf['datasets']['sidis']['filters']={}
conf['datasets']['sidis']['filters'][0]={}
#conf['datasets']['sidis']['filters'][0]['idx']=[1000,1001,1004,1005]
conf['datasets']['sidis']['filters'][0]['idx']=[1000,1001,1004,1005,1002,1003,1006,1007]
#conf['datasets']['sidis']['filters'][0]['idx']=[1000,1001,1004,1005,1002,1003,1006,1007,5001,5002]
#conf['datasets']['sidis']['filters'][0]['idx']=[5001,5002]
#conf['datasets']['sidis']['filters'][0]['filter']="z>0.2 and z<0.6 and Q2>1.69 and (pT/z)**2<0.25*Q2 and 2*yh<-0.69" # R filter exp(2 yh)<0.5 npts    = 368 chi2    = 589.505730
#conf['datasets']['sidis']['filters'][0]['filter']="z>0.2 and z<0.6 and Q2>1.69 and  pT>0.2 and pT<0.9" # no R filter npts    = 807 chi2    = 932.408505 
#conf['datasets']['sidis']['filters'][0]['filter']="z>0.2 and z<0.6 and Q2>1.69 and (pT/z)**2<0.25*Q2 and yp-yh>2." # rapidity difference npts    = 481 chi2 = 683.694846    = 683.694846
conf['datasets']['sidis']['filters'][0]['filter']="z>0.2 and z<0.6 and Q2>1.69 and (pT/z)**2<0.15*Q2 and yp-yh>3." # rapidity difference npts    = 481 chi2 = 683.694846    = 683.694846
#conf['datasets']['sidis']['filters'][0]['filter']="z<0.6 and Q2>1.69 and (pT/z)**2<0.25*Q2 and yp-yh>0." #npts    = 513 chi2    = 722.943137
#conf['datasets']['sidis']['filters'][0]['filter']="z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9" # npts    = 978 chi2    = 1145.231816

conf['datasets']['sidis']['xlsx']={}

conf['datasets']['sidis']['xlsx'][1000]='../database/sidis/expdata/1000.xlsx'  # |  proton   | pi+    | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1001]='../database/sidis/expdata/1001.xlsx'  # |  proton   | pi-    | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1004]='../database/sidis/expdata/1004.xlsx'  # |  deuteron | pi+    | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1005]='../database/sidis/expdata/1005.xlsx'  # |  deuteron | pi-    | M_Hermes | hermes
#
conf['datasets']['sidis']['xlsx'][1002]='../database/sidis/expdata/1002.xlsx'  # |  proton   | k+    | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1003]='../database/sidis/expdata/1003.xlsx'  # |  proton   | k-    | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1006]='../database/sidis/expdata/1006.xlsx'  # |  deuteron | k+    | M_Hermes | hermes
conf['datasets']['sidis']['xlsx'][1007]='../database/sidis/expdata/1007.xlsx'  # |  deuteron | k-    | M_Hermes | hermes

#conf['datasets']['sidis']['xlsx'][5001]='../database/sidis/expdata/5001.xlsx'  # |  deuteron | h+    | M_Compass| compass
#conf['datasets']['sidis']['xlsx'][5002]='../database/sidis/expdata/5002.xlsx'  # |  deuteron | h-    | M_Compass| compass


conf['datasets']['sidis']['norm']={}
for k in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1}






