
conf = {}

def load_config(fname):
    L=open(fname).readlines()
    D = {}
    for l in L: exec l.replace('<<','').replace('>>','') in D
    conf.update(D['conf'])
