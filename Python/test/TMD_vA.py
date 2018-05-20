#!/usr/bin/env python
import sys, os
import numpy as np
from scipy.special import j0, j1
import scipy.integrate as integrate
#import fCJLO, fCJNLO
from external.CJLIB.CJ import CJ

def alphas(Q):
    mb = 4.18
    mc = 1.28
    if (Q < mc):
        nf = 3
        lam = 0.372318
    elif (Q < mb):
        nf = 4
        lam = 0.326514
    else:
        nf = 5
        lam = 0.22878
    b0 = 11.0 - 2.0 / 3.0 * nf
    b1 = 102.0 - 38.0 / 3.0 * nf
    tt = 2.0 * np.log(Q/lam)
    result = 4.0 * np.pi / (b0 * tt) * (1.0 - b1 * np.log(tt) / (b0 * b0 * tt))
    return result

def bessel0(z):
    return j0(z)

def coll_int(z, zv, Q, ip):
    alfpi = alphas(Q) / np.pi
    Pqq = alfpi * CF * (1.0 - z) / (2.0 * z)
    Pgq = alfpi * CF * (z / 2.0 + (1.0 + (1.0 - z) / 2.0) / z * np.log(z)) / z
    Pqg = alfpi * TR * (z * (1.0 - z)) / z
    Pgg = alfpi * 2.0 * CA * (z / (1.0 - z) + (1.0 - z) / z + z * (1.0 - z)) * np.log(z) / z
    f = CT14Demo(zv / z, Q)
    if ip == 0:
        qsum = f[1] + f[-1] + f[2] + f[-2] + f[3] + f[-3] + f[4] + f[-4] + f[5] + f[-5]
        res = Pgg * f[0] + Pqg * qsum
    else:
        res = Pqq * f[ip] + Pgq * f[0]
    return res

def coll_ff(iq, z, Q):
    f = CT14Demo(z, Q)
    fqgcommon, err = integrate.quad(coll_int, z, 1, args = (z, Q, iq))
    return f[iq] + fqgcommon



def f_tilde(b, x, Q, iq, order):
    c0 = 1.22919
    Q0 = c0 / b
    Q0min = c0 / bmax
    if Q0 < 0.5:
        Q0 = 0.5
    if iq == 0:
        FNP = 1.0
        res = coll_ff(0, x, Q0) * kernel_g(b, Q0, Q, Q0, Q, order) * FNP
    elif b <= bmax:
        FNP = 1
        res = coll_ff(iq, x, Q0) * kernel_q(b, Q0, Q, Q0, Q, order) * FNP
    elif b > bmax:
        h = 1e-6
        fx = f_tilde(bmax, x, Q, iq, order)
        fxh = f_tilde(bmax - h, x, Q, iq, order)
        fx2h = f_tilde(bmax - 2.0 * h, x, Q, iq, order)
        g2 = 0.9
        db1 = (fx - fxh) / h
        db2 = (fx2h - 2.0 * fxh + fx) / (h ** 2)
        alp = 0.5 * (1.0  / (db1 / (2.0 * fx * bmax) + g2) * (g2 + db2 / (2.0 * fx) - db1 ** 2 / (2.0 * fx ** 2)) + 1)
        g1 = -bmax ** (1.0 - 2.0 * alp) / alp * (db1 / (2.0 * fx) + g2 * bmax)
        FNP = np.exp(-g1 * (b ** (2.0 * alp) - bmax ** (2.0 * alp)) - g2 * (b ** 2 - bmax ** 2))
        res = coll_ff(iq, x, Q0min) * kernel_q(iq, x, Q0min, Q, Q0min, Q, order) * FNP
    else:
        return 0
    return res

def tmd_int(b, kt, x, Q, iq, order):
    return f_tilde(b, x, Q, iq, order) * b / (2.0 * np.pi) * bessel0(b * kt)

def tmd(kt, x, Q, iq, order):
    res, err = integrate.quad(tmd_int, 0, np.inf, args = (kt, x, Q, iq, order))
    return res
        
    
        









if __name__=='__main__':

    conf={}
    conf['path2CJ']='./'
    conf['order']='NLO'
    cjNLO=CJ(conf)

    print cjNLO.get_f(0.1, 1.62)
    sys.exit()
    
