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
conf['params']['pdf']['widths0 valence']  = {'value':<<    5.50013020807587293959e-01>>,'fixed':False,'min':0,'max':10}
conf['params']['pdf']['widths0 sea']      = {'value':<<    1.06392801081480636860e-01>>,'fixed':False,'min':0,'max':10}

conf['params']['ff']={}
conf['params']['ff']['widths0 pi+ fav']   = {'value':<<    1.18205757936943378628e-01>>,'fixed':False,'min':0,'max':1}
conf['params']['ff']['widths0 pi+ unfav'] = {'value':<<    5.64424548261521108472e-01>>,'fixed':False,'min':0,'max':1}
conf['params']['ff']['widths0 k+ fav']    = {'value':<<    1.33775551334897768729e-01>>,'fixed':False,'min':0,'max':1}
conf['params']['ff']['widths0 k+ unfav']  = {'value':<<    1.83430383139499497691e-01>>,'fixed':False,'min':0,'max':1}

conf['params']['gk']={}
conf['params']['gk']['g2']   = {'value':<<    1.31574342618341721955e-01>>,'fixed':False,'min':0,'max':10}
conf['params']['gk']['Q0'] = {'value':<<    1.00000000000000000000e+00>>,'fixed':True,'min':1,'max':10}
conf['params']['gk']['bmax'] = {'value':<<    1.00000000000000000000e+00>>,'fixed':True,'min':0.01,'max':10}


############################################################################
# set data sets

conf['datasets']={}

# SIDIS

conf['datasets']['sidis']={}


conf['datasets']['sidis']['filters']={}
conf['datasets']['sidis']['filters'][0]={}
#conf['datasets']['sidis']['filters'][0]['idx']=[1000,1001,1004,1005,1002,1003,1006,1007]
#conf['datasets']['sidis']['filters'][0]['idx']=[1000,1001,1004,1005,1002,1003,1006,1007,5001,5002]
conf['datasets']['sidis']['filters'][0]['idx']=[5001,5002]
#conf['datasets']['sidis']['filters'][0]['filter']="z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9"
#conf['datasets']['sidis']['filters'][0]['filter']="z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9 and (pT/z)**2<0.5*Q2 and yp-yh>2."
conf['datasets']['sidis']['filters'][0]['filter']="z<0.9 and Q2>1.69 and (pT/z)**2<0.25*Q2 and yp-yh>2."

conf['datasets']['sidis']['xlsx']={}

#conf['datasets']['sidis']['xlsx'][1000]='../database/sidis/expdata/1000.xlsx'  # |  proton   | pi+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1001]='../database/sidis/expdata/1001.xlsx'  # |  proton   | pi-    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1004]='../database/sidis/expdata/1004.xlsx'  # |  deuteron | pi+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1005]='../database/sidis/expdata/1005.xlsx'  # |  deuteron | pi-    | M_Hermes | hermes 
#
#conf['datasets']['sidis']['xlsx'][1002]='../database/sidis/expdata/1002.xlsx'  # |  proton   | k+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1003]='../database/sidis/expdata/1003.xlsx'  # |  proton   | k-    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1006]='../database/sidis/expdata/1006.xlsx'  # |  deuteron | k+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1007]='../database/sidis/expdata/1007.xlsx'  # |  deuteron | k-    | M_Hermes | hermes 

conf['datasets']['sidis']['xlsx'][5001]='../database/sidis/expdata/5001.xlsx'  # |  deuteron | h+    | M_Compass| compass 
conf['datasets']['sidis']['xlsx'][5002]='../database/sidis/expdata/5002.xlsx'  # |  deuteron | h-    | M_Compass| compass 
    
    
conf['datasets']['sidis']['norm']={}
for k in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1} 






