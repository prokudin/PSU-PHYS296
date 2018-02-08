Adaptive Ogata algorithm. To run use the command

from tools.hankel.inverters import Ogata
ogata = Ogata()
ogata.adog(w, qT, nu, Nmax, bmin, bmax)

-w is a function of some variable bT to be inverted. 
-qT is the Fourier conjugate variable to bT. 
-nu is typical integer specifying the bessel function. 
-Nmax is the maximum number of calls to the function on the LAST interation. Total number of function calls will go like ~10*Nmax
-bmin and bmax are used to find the maximum of w function in b space. The program will call the function 20 times in this region seeking the maximum value of w(b) and the bcrit that corresponds to this. Specify bmin just below bcrit and bmax just above bcrit. 
