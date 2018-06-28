import numpy as np
from numpy import log

####################################################
##########Interacciones tipo P2P####################
####################################################

def P2P(m,y,x):
    return -sum(m[:]*log(abs(x-y[:])))
    #return sum(m[:]/abs(x-y[:]))


def P2M(m,fact,yp,ym):    
    p = len(fact)
    #np = len(yp)    
    M = np.zeros(p, dtype=complex)
    for n in range(p):
        M[n] = sum(((ym-yp[:])**n)*m[:]/fact[n])
    #print M
    return M

def M2M(Mc, fact, yM, yc):
    p = len(Mc)
    M = np.zeros(p, dtype=complex)
    for n in range(p):
        for k in range(n+1):
            M[n] += Mc[k]*((yM-yc)**(n-k))/fact[n-k]
    return M
    
def M2L(M,fact,yM,yL):
    p = len(M)
    L = np.zeros(p, dtype=complex)
    for k in range(p):
        for n in range(p-k):
            #L[k] += ((-1.)**(n+k))*fact[n+k]*((yL-yM)**(n+k))*(abs(yL-yM)**(-2*n-2*k-1))*M[n]
            if n+k == 0:
                L[k] += -log(abs(yL-yM))*M[n]
            else:
                L[k] += (((-1.)**(n+k))*(fact[n+k-1])*(yL-yM)**(-n-k))*M[n]
        #print L[k]
    return L

def L2L(L,fact,xL,xc):
    p = len(L)
    Lc = np.zeros(p, dtype=complex)
    for n in range(p):
        for k in range(n,p):
            Lc[n] += ((xc-xL)**(k-n))*L[k]/fact[k-n]
    return Lc
            
def L2P(L,fact,xp,xL):
    p = len(L)
    phi = 0.0    
    for k in range(p):
        phi += ((xp-xL)**k)*L[k]/fact[k]
        #print phi
    return phi
