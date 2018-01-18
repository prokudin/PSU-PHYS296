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


conf['params']['ppdf']={}
conf['params']['ppdf']['widths0 valence']  = {'value':<<    1.09361479714813469855e-01>>,'fixed':False,'min':0,'max':10}
conf['params']['ppdf']['widths0 sea']      = {'value':<<    2.18478703409406538327e+00>>,'fixed':False,'min':0,'max':10}


conf['params']['pdf']={}
conf['params']['pdf']['widths0 valence']  = {'value':<<    2.50000000000000000000e-01>>,'fixed':True,'min':0,'max':1}
conf['params']['pdf']['widths0 sea']      = {'value':<<    2.50000000000000000000e-01>>,'fixed':True,'min':0,'max':1}

conf['params']['ff']={}
conf['params']['ff']['widths0 pi+ fav']   = {'value':<<    2.00000000000000011102e-01>>,'fixed':True,'min':0,'max':1}
conf['params']['ff']['widths0 pi+ unfav'] = {'value':<<    2.00000000000000011102e-01>>,'fixed':True,'min':0,'max':1}
conf['params']['ff']['widths0 k+ fav']    = {'value':<<    2.00000000000000011102e-01>>,'fixed':True,'min':0,'max':1}
conf['params']['ff']['widths0 k+ unfav']  = {'value':<<    2.00000000000000011102e-01>>,'fixed':True,'min':0,'max':1}


conf['params']['gk']={}
conf['params']['gk']['g2']   = {'value':<<    0.00000000000000000000e+00>>,'fixed':True,'min':0,'max':10}
conf['params']['gk']['Q0'] = {'value':<<    1.67398382026420211588e+00>>,'fixed':True,'min':1,'max':10}

############################################################################
# set data sets

conf['datasets']={}

# SIDIS

conf['datasets']['sidis']={}


conf['datasets']['sidis']['filters']={}
conf['datasets']['sidis']['filters'][0]={}
conf['datasets']['sidis']['filters'][0]['idx']=[2000,2001,2003,2004]
#conf['datasets']['sidis']['filters'][0]['idx']=[2000,2001]
#conf['datasets']['sidis']['filters'][0]['filter']="z<0.6 and Q2>1.69 and pT>0.2 and pT<0.9"
conf['datasets']['sidis']['filters'][0]['filter']="z>0"

conf['datasets']['sidis']['xlsx']={}

#conf['datasets']['sidis']['xlsx'][1000]='../database/sidis/expdata/1000.xlsx'  # |  proton   | pi+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1001]='../database/sidis/expdata/1001.xlsx'  # |  proton   | pi-    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1004]='../database/sidis/expdata/1004.xlsx'  # |  deuteron | pi+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1005]='../database/sidis/expdata/1005.xlsx'  # |  deuteron | pi-    | M_Hermes | hermes 

#conf['datasets']['sidis']['xlsx'][1002]='../database/sidis/expdata/1002.xlsx'  # |  proton   | k+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1003]='../database/sidis/expdata/1003.xlsx'  # |  proton   | k-    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1006]='../database/sidis/expdata/1006.xlsx'  # |  deuteron | k+    | M_Hermes | hermes 
#conf['datasets']['sidis']['xlsx'][1007]='../database/sidis/expdata/1007.xlsx'  # |  deuteron | k-    | M_Hermes | hermes 

    
conf['datasets']['sidis']['xlsx'][2000]='../database/sidis/expdata/2000.xlsx'  # |  proton | pi+    | ALL | clas 
conf['datasets']['sidis']['xlsx'][2001]='../database/sidis/expdata/2001.xlsx'  # |  proton | pi-    | ALL | clas 
#conf['datasets']['sidis']['xlsx'][2002]='../database/sidis/expdata/2002.xlsx'  # |  proton | pi0    | ALL | clas 
    
conf['datasets']['sidis']['xlsx'][2003]='../database/sidis/expdata/2003.xlsx'  # |  proton | pi+    | ALL | compass 
conf['datasets']['sidis']['xlsx'][2004]='../database/sidis/expdata/2004.xlsx'  # |  proton | pi-    | ALL | compass 
    
conf['datasets']['sidis']['norm']={}
for k in conf['datasets']['sidis']['xlsx']: conf['datasets']['sidis']['norm'][k]={'value':1,'fixed':True,'min':0,'max':1} 






