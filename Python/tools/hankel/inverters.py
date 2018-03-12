#!/usr/bin/env python2

import sys,os
import numpy as np
import heapq
#import vegas
from scipy.special import jv, jn_zeros, yv
from scipy.optimize import fsolve
from scipy.integrate import quad, fixed_quad
import warnings
warnings.filterwarnings('ignore')

class Ogata:

    def ogata(self, f,h,N, nu):
        zeros=jn_zeros(nu,N)
        xi=zeros/np.pi
        Jp1=jv(nu+1,np.pi*xi)
        w=yv(nu, np.pi * xi) / Jp1
        get_psi=lambda t: t*np.tanh(np.pi/2*np.sinh(t))
        get_psip=lambda t:np.pi*t*(-np.tanh(np.pi*np.sinh(t)/2)**2 + 1)*np.cosh(t)/2 + np.tanh(np.pi*np.sinh(t)/2)
        knots=np.pi/h*get_psi(h*xi)
        Jnu=jv(nu,knots)
        psip=get_psip(h*xi)
        F=f(knots)
        return np.pi*np.sum(w*F*Jnu*psip)
    
    def find_bmax(self, f, bmin, bmax, bnum = 10, iterations = 2, first = True):
        _iterations = iterations-1
        if first == True:
            bT = np.exp(np.linspace(np.log(bmin), np.log(bmax), bnum))
        else:
            bT = np.linspace(bmin, bmax, bnum)
        f_bT = [f(b) for b in bT]
        fb2max = heapq.nlargest(2, f_bT)[1]
        imax = [i for i, j in enumerate(f_bT) if j == max(f_bT)][0]
        i2max = [i for i, j in enumerate(f_bT) if j == fb2max][0]
        b_list_max = bT[imax]
        b_list_2_max = bT[i2max]
        #print b_list_max
        if _iterations >= 0:
            return self.find_bmax(f, b_list_2_max, b_list_max, bnum = 10, iterations = _iterations, first = False)
        else:
            return b_list_max
    
#    def get_ogata_params(self, w, xmin, xmax, qT, nu):
#        zero1 = jn_zeros(nu, 1)[0]
#        h = fsolve(lambda h: xmin-zero1*np.tanh(np.pi/2*np.sinh(h/np.pi*zero1)), xmin)[0]
#        k = fsolve(lambda k: xmax-np.pi*k*np.tanh(np.pi/2*np.sinh(h*k)), xmax)[0]
#        N = int(k)
#        return h,k,N

#    def Ogatainversion(self, w, xmin, xmax, qT, nu):
#        h,k,N=self.get_ogata_params(w, xmin, xmax, qT, nu)
#        return 1/(2*np.pi)*self.ogata(lambda x: w(x/qT)/qT,h,N, nu)
    
    def get_ogata_params_b(self, w, bmin, bmax, qT, nu):
        zero1 = jn_zeros(nu, 1)[0]
        h = fsolve(lambda h: bmin-zero1/qT*np.tanh(np.pi/2*np.sinh(h/np.pi*zero1)), bmin)[0]
        k = fsolve(lambda k: bmax-np.pi/qT*k*np.tanh(np.pi/2*np.sinh(h*k)), bmax)[0]
        if k<0:
            k = fsolve(lambda k: bmax-np.pi/qT*k*np.tanh(np.pi/2*np.sinh(h*k)), -k)[0]
        N = int(k+1)
        return h,N
    
    def get_b_vals(self, bc, ib, it):
        ccut = 2.0
        #Original b range
        bmin = bc/((ccut)**ib)
        bmax = bc*ccut**it
        #Double bmax
        bmin1 = bc/((ccut)**ib)
        bmax1 = bc*ccut**it*ccut
        #Half bmin
        bmin2 = bc/((ccut)**ib)/ccut
        bmax2 = bc*ccut**it
        #Half bmin and double bmax
        bmin3 = bc/((ccut)**ib)/ccut
        bmax3 = bc*ccut**it*ccut
        return bmin, bmax, bmin1, bmax1, bmin2, bmax2, bmin3, bmax3
    
    def compare(self, w, bmin, bmax, bmin1, bmax1, bmin2, bmax2, bmin3, bmax3, qT, nu, lib):
        #Original b range
        h,N=self.get_ogata_params_b(w, bmin, bmax, qT, nu)
        #Double bmax
        h1,N1=self.get_ogata_params_b(w, bmin1, bmax1, qT, nu)
        #Half bmin
        h2,N2=self.get_ogata_params_b(w, bmin2, bmax2, qT, nu)
        #Half bmin and double bmax
        h3,N3=self.get_ogata_params_b(w, bmin3, bmax3, qT, nu)
        try:
            val1 = lib[str(bmin1)+','+str(bmax1)]
        except:
            val1 = 1/(2*np.pi)*self.ogata(lambda x: w(x/qT)/qT,h1,N1, nu)
            lib[str(bmin1)+','+str(bmax1)] = val1
        try:
            val2 = lib[str(bmin2)+','+str(bmax2)]
        except:
            val2 = 1/(2*np.pi)*self.ogata(lambda x: w(x/qT)/qT,h2,N2, nu)
            lib[str(bmin2)+','+str(bmax2)] = val2
        try:
            val3 = lib[str(bmin3)+','+str(bmax3)]
        except:
            val3 = 1/(2*np.pi)*self.ogata(lambda x: w(x/qT)/qT,h3,N3, nu)
            lib[str(bmin3)+','+str(bmax3)] = val3
        return abs((val1-val3)/val3)<abs((val2-val3)/val3), lib
    
    def adog(self, w, qT, nu, Nmax, bmin, bmax, itmax=100, ib = 1, it = 1, bpeak = [], lib = None, bminlast = 0, bmaxlast = 0):
        if lib is None:
            lib = {}
        if ib+it< itmax:
            _lib = lib
            if not bpeak:
                bc = self.find_bmax(w, bmin, bmax)
            else:
                bc = bpeak[0]
            bmin, bmax, bmin1, bmax1, bmin2, bmax2, bmin3, bmax3 = self.get_b_vals(bc, ib, it)
            h2, N2 = self.get_ogata_params_b(w, bmin2, bmax2, qT, nu)
            if N2<Nmax:
                _bool, __lib = self.compare(w, bmin, bmax, bmin1, bmax1, bmin2, bmax2, bmin3, bmax3, qT, nu, _lib)
                if _bool == True:
                    return self.adog(w, qT, nu, Nmax, bmin, bmax, itmax = 100, ib = ib+0, it = it+1, bpeak = [bc], lib = __lib, bminlast = bmin3, bmaxlast = bmax3)
                else:
                    return self.adog(w, qT, nu, Nmax, bmin, bmax, itmax = 100, ib = ib+1, it = it+0, bpeak = [bc], lib = __lib, bminlast = bmin3, bmaxlast = bmax3)
            else:
                return _lib[str(bminlast)+','+str(bmaxlast)]
        else:
            print('itmax reached')
    def adog2(self, w, qT, nu, Nmax, Q, itmax=100, ib = 1, it = 1, bpeak = [], lib = None, bminlast = 0, bmaxlast = 0):
        if lib is None:
            lib = {}
        if ib+it< itmax:
            _lib = lib
            bc = 1/Q
            bmin, bmax, bmin1, bmax1, bmin2, bmax2, bmin3, bmax3 = self.get_b_vals(bc, ib, it)
            h2, N2 = self.get_ogata_params_b(w, bmin2, bmax2, qT, nu)
            if N2<Nmax:
                _bool, __lib = self.compare(w, bmin, bmax, bmin1, bmax1, bmin2, bmax2, bmin3, bmax3, qT, nu, _lib)
                if _bool == True:
                    return self.adog2(w, qT, nu, Nmax, Q, itmax = 100, ib = ib+0, it = it+1, bpeak = [bc], lib = __lib, bminlast = bmin3, bmaxlast = bmax3)
                else:
                    return self.adog2(w, qT, nu, Nmax, Q, itmax = 100, ib = ib+1, it = it+0, bpeak = [bc], lib = __lib, bminlast = bmin3, bmaxlast = bmax3)
            else:
                return _lib[str(bminlast)+','+str(bmaxlast)]
        else:
            print('itmax reached')
class Quad:

    def quadinv(self, w, q, nu, eps):
        quadreturn = quad(lambda bT: jv(nu,q*bT)*w(bT),0.0,5.0, epsabs = 0.0, epsrel = eps)
        return 1/(2*np.pi)*quadreturn[0], 1/(2*np.pi)*quadreturn[1]

class Fix_Quad:

    def fix_quadinv(self, w, q, nu, num):
        quadreturn = fixed_quad(lambda bT: jv(nu,q*bT)*w(bT),0.0,2.0, n = num)
        return 1/(2*np.pi)*quadreturn[0]#, 1/(2*np.pi)*quadreturn[1]

#class Vegas:
#
#    def transform(self, f, p):
#        return f(np.tan(p))*(1/np.cos(p))**2
#
#    def MCinv(self, f, q, nu, m):
#        integ = vegas.Integrator([[0, np.pi/2.0]])
#        result = integ(lambda p: 1/(2*np.pi)*self.transform(f, p)*jv(nu, q*np.tan(p)), nitn=10, neval=int(m))[0]
#        lst = str(result).replace('+-', '(').replace(')', '(').split('(')
#        num_zeros = 0
#        if str(lst[0][0]) == '-':
#            num_zeros = -1
#        if len(lst[1]) == 3:
#            integral = float(lst[0])
#            error = float(lst[1])
#        else:
#            num_zeros = num_zeros+len(lst[0])-4
#            lst[1] = '0.'+lst[1].zfill(num_zeros+2)
#            integral = float(lst[0])
#            error = float(lst[1])
#        return integral, error



if __name__== "__main__":
    ogata = Ogata()
    adquad = Quad()
#    vegasmc = Vegas()
    print ogata.adog(lambda x: x*np.exp(-x**2), 1.0, 0, 100, 0.5, 1.0)
    print adquad.quadinv(lambda x: x*np.exp(-x**2), 1.0, 0)
#    print vegasmc.MCinv(lambda x: x*np.exp(-x**2), 1.0, 0, 1000)
