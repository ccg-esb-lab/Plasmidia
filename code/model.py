import numpy as np
import scipy as scipy
import scipy.linalg
import math
from decimal import Decimal
from decimal import *


from params import *


def alpha_star(n,thiscc=pBGT_cost):
    astar=Decimal(thiscc*n)+mutationrate(n)*Decimal((1-thiscc*n))
    return astar

def mutationrate(n):
    #print(r)
    ro = ρN*n
    rod = Decimal(ρN*n)
    rd=Decimal(r)
    σd=Decimal(σ)
    aux=Decimal(2**-n)
    ui=1-( ((rd*aux)+rod)/( (rd*aux)*np.exp((rd*aux+rod)*σd) +rod) )
    return ui

 #iterates t days the non-antibiotics equation
def solucionvectorial(x0, n, c, t,u):
    X=[None]*(t+1)
    X[0] = x0
    for i in range(t):
        x0=X[i]
        X[i+1] = (1-u)*(1-c)*x0/( ((1-c)*x0)+1-x0)
    return X[1:]

# non-antibiotics equation
def solucionAntibiotic(x0, n, c, a,u):
    #print("a",a,x0,"x")
    if (a==1. and x0==0.):
        y=np.nan
    else:
        y=(u-1)*(c-1)*x0/(1-a+x0*(a-c))
    return y

#Alternates T-1 days non-antibiotics with 1 antiobiotics
def solutionPlasmidsPeaks(x0, n, c, Tmax, T, a):
    u = mutationrate(n)
    c=Decimal(c)
    a=Decimal(a)
    x0=Decimal(x0)
    
    X=[None]*(Tmax+1)
    X[0] = x0
    if (T==0):
        p=0
    else:
        p = Tmax/T #number of peaks
    for k in range(int(p)):
        X[k*T+1:(k+1)*T] = solucionvectorial(X[k*T], n, c, T-1,u)
        X[(k+1)*T] = solucionAntibiotic(X[(k+1)*T-1] , n, c, a,u)
    X[int(p)*T+1:] = solucionvectorial(X[int(p)*T], n, c, len(X[int(p)*T+1:]),u)
    return X





######  additionals
def get_first_rep(iX):
    iX=[float(x) for x in iX]
    indxxs=[]
    l=len(iX)
    e=int(l/10)
    iX2=iX[e:]
    if(e>=10000):
        indxxs=list(np.where(np.isin(iX,iX[-1]))[0])
        if(indxxs==[l-1]):
            indxxs=[]
    else:
        for xi,x in enumerate(iX2): 
            if ((x in iX2[xi+1:])): 
            #indxxs.append(xi)
                indxxs=[xi+e]
                break 
    
    return indxxs


def solutionPlasmidsOnePeak(x0, n, c, Tmax, T, a):
    u = mutationrate(n)
    c=Decimal(c)
    a=Decimal(a)
    x0=Decimal(x0)
    X=[None]*(Tmax+1)
    X[0] = x0
    if (T==0):
        p=0
    else:
        p = Tmax/T #number of peaks
    
    
    f=True
    for k in range(int(p)):
        #print(kk,end=", ")
        
        X[k*T+1:(k+1)*T] = solucionvectorial(X[k*T], n, c, T-1,u)
        if(f):
            X[(k+1)*T] = solucionAntibiotic(X[(k+1)*T-1] , n, c, a,u)
            f=False
        else:
            X[(k+1)*T] = solucionAntibiotic(X[(k+1)*T-1] , n, c, 0,u)
            #X[k*T+1:(k+1)*T] = solucionvectorial(X[k*T], n, c, T-1,u)
    X[int(p)*T+1:] = solucionvectorial(X[int(p)*T], n, c, len(X[int(p)*T+1:]),u)
    return X


    
def get_first_rep_peak(iX,tT,tnts):
    iX=[float(x) for x in iX]
    indxxs=[]
    iX2=[]#iX[e:]
    for ti in range(0,tnts+1):
        iX2.append(iX[ti*tT])
    #print(iX)
    #print(iX2)
    l=len(iX2)
    
    if(l>=10000):
        indxxs=list(np.where(np.isin(iX2,iX2[-1]))[0])
        indxxs=[x*tT for x in indxxs]
        if(indxxs==[l-1]):
            indxxs=[]
        
    else:
        for xi,x in enumerate(iX2): 
            if ((x in iX2[xi+1:])): 
                indxxs=[xi*tT]
                break 
    #print(indxxs)
    return indxxs

def get_fancy_map(x00,iX,tnts,tT):
    sX=[x00,x00]
    sX2=[x00,iX[1]]
    for i in range(1,tnts+1):
        vi=i*tT
#        vi2=i*tT
        
        sX.append(iX[vi-2])
        sX.append(iX[vi-1])
        sX.append(iX[vi])
        sX2.append(iX[vi-1])
        sX2.append(iX[vi])
        sX2.append(iX[vi+1])    
    return sX,sX2

def get_fancy_map_avg(x00,iX,tnts,tT):
    v1=(iX[tT]+iX[tT-1])/2
    sX=[x00,x00]
    sX2=[x00,v1]
    sX=[x00]
    sX2=[x00]
    for i in range(1,tnts):
        vi=i*tT
        v1=(iX[vi]+iX[vi-1])/2
        #print((vi,iX[vi],iX[vi-1],v1),end=",")
        vi2=(i+1)*tT
        v2=(iX[vi2]+iX[vi2-1])/2
        sX.append(v1)
        sX2.append(v2)
        
    return sX,sX2



def f0(x,g,mu):
    y=x*(1-g)/(x*(1-g)+(1-x))*(1-mu)
    return y

def f1(x,gamma,alpha,mu):
    y=x*((1-gamma)/((1-alpha)+x*(alpha-gamma)))*(1-mu)
    return y
epsg=1e-4


