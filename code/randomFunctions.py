import sys
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as scipy
import scipy.linalg
from itertools import *
import re
from random import *
import seaborn as sns
from shapely import *
from shapely.geometry import *
from descartes.patch import PolygonPatch
import math,string,random
import matplotlib.patches as mpatches
from urllib.request import urlopen
from scipy.stats import binom
import json
import pylab as P
from shapely import *
from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import *
from descartes import PolygonPatch

from model import *

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def jsonSaveEs(fnamej,tosaveEs):
    
    with open(fnamej, mode='w') as file:
        json.dump(tosaveEs,file,cls=NpEncoder)
    file.close()
    
def jsonLoadEs(fnamej):
    with open(fnamej, 'r') as file:
        ESj = json.load(file)
    file.close()
    return ESj


def invert_Ev(Ev):
    Ev=np.where(Ev==1,-1,Ev)
    Ev=np.where(Ev==0, 1, Ev)
    Ev=np.where(Ev==-1, 0, Ev)
    return Ev


def solutionPlasmidsEnv(x0, n, c, Tmax, E, a):
    u = mutationrate(n)
    c=Decimal(c)
    a=Decimal(a)
    x0=Decimal(x0)
    #c=round(c,5)
    X=[None]*(Tmax+1)
    X[0] = x0
    
    for ti,Ta in enumerate(E):
        if(Ta==0):
            X[ti+1] = solucionAntibiotic(X[ti] , n, c, 0,u)
        else:
            X[ti+1] = solucionAntibiotic(X[ti] , n, c, a,u)
    
    return X




def H2(data):
    tdict={}
    entropy = 0
    for x in data:
        tdict[x]=0
    for x in data:
        tdict[x]+=1
    #print(tdict)
    for k,v in tdict.items():
        p_x = v/(len(data))
        if p_x > 0:
            entropy += -p_x*math.log(p_x, 2)
    #entropy=-entropy
    return entropy

def getShannon(Es):
    shannon=[]
    for E in Es:
        shannon.append(E['shannon'])
    return shannon

def getFraction(Es):
    ff=[]
    ln=len(Es[0]['Ev'])
    for E in Es:
        f=sum(E['Ev'])/ln
        ff.append(f)
    return ff

def getRandomEnvironment(T=100, p=0.5, alphabet=[0,1]):
    n=max(alphabet)
    EV= binom.rvs(n, p, size=T)
    E=dict()
    E['Ev']=[x for x in EV]
    E['shannon']=H2(E['Ev'])

    return E

def plotEnvironmentStats_Bins(Es,numBins,Es_BinsH,Es_BinsF):
    f, ax = plt.subplots(1,4,figsize=(16,4))
    f.set_facecolor('white')
    Est=[]
    for ib,this_Es_byBin in enumerate(Es_BinsH):
        for ie,this_Eindex in enumerate(this_Es_byBin):
            thisE=Es[this_Eindex]
            Est.append(thisE)
    #Plots shannon distribution H
    n3, bins3, patches3 = ax[0].hist(getShannon(Est), bins=numBins )
    P.setp(patches3, 'facecolor', 'green', 'alpha', 0.75)
    ax[0].set_xlabel('Shannon Entropy')
    ax[0].set_ylabel('')

    #Plots Fraction distribution distribution H
    n3, bins3, patches3 = ax[1].hist(getFraction(Est), bins=numBins )
    P.setp(patches3, 'facecolor', 'purple', 'alpha', 0.75)
    ax[1].set_xlabel('Fraction of 1s')
    ax[1].set_ylabel('')
    
    Est=[]
    for ib,this_Es_byBin in enumerate(Es_BinsF):
        for ie,this_Eindex in enumerate(this_Es_byBin):
            thisE=Es[this_Eindex]
            Est.append(thisE)

    #Plots shannon distribution F
    n3, bins3, patches3 = ax[2].hist(getShannon(Est), bins=numBins )
    P.setp(patches3, 'facecolor', 'green', 'alpha', 0.75)
    ax[2].set_xlabel('Shannon Entropy')
    ax[2].set_ylabel('')

    #Plots Fraction distribution distribution F
    n3, bins3, patches3 = ax[3].hist(getFraction(Est), bins=numBins )
    P.setp(patches3, 'facecolor', 'purple', 'alpha', 0.75)
    ax[3].set_xlabel('Fraction of 1s')
    ax[3].set_ylabel('')

def plotEnvironmentStats(Est,numBins):
    f, ax = plt.subplots(1,2,figsize=(12,4))
    f.set_facecolor('white')
   #Plots shannon distribution
    n3, bins3, patches3 = ax[0].hist(getShannon(Est), bins=numBins )
    P.setp(patches3, 'facecolor', 'green', 'alpha', 0.75)
    ax[0].set_xlabel('Shannon Entropy')
    ax[0].set_ylabel('')

    #Plots Fraction distribution distribution
    n3, bins3, patches3 = ax[1].hist(getFraction(Est), bins=numBins )
    P.setp(patches3, 'facecolor', 'purple', 'alpha', 0.75)
    ax[1].set_xlabel('Fraction of 1s')
    ax[1].set_ylabel('')

def getRandomEnvironmentsUniformEntropies(ndays,max_tries,numBins,numEs):
    
    H_bins=np.arange(0,1+1/numBins,1/numBins)
    print("Bins are : ",H_bins)
    min_H=(-1/ndays*math.log(1/ndays, 2))+(-(ndays-1)/ndays*math.log((ndays-1)/ndays, 2))
    print("Min H :",min_H)
    max_frac_completeH=1
    if(min_H>H_bins[1]):
        print("The bin size is smaller than the lowest possible Entropy")
        print("Some bins will not be covered")
        print("Do you wish to continue?")
        pchoice = input("Yes(y)/No(n)")
        if(pchoice=="n"):
            return [],[],[]
        elif(pchoice=="y"):
            print("OK!")
            max_frac_completeH=1-1/numBins
            max_frac_complete=1-1/numBins
        else:
            print("All you base are belong to us")
            return 0;



    H_bins_countF=[int(x) for x in np.zeros(len(H_bins))]
    H_bins_countH=[int(x) for x in np.zeros(len(H_bins))]
    nbins=len(H_bins)
    nbins=2
    Es=[]
    Es_byBins0H=[[] for x in H_bins]
    Es_byBins0F=[[] for x in H_bins]

    H_bins_completedH=np.zeros(len(H_bins))
    H_bins_completedF=np.zeros(len(H_bins))
    
    frac_completeH=0
    frac_completeF=0
    tshannon=0
    env_index=0
    cycle=1
    fA=0
    rho=1
    
    last_sum=0
    cycle_buffer=1
    while((cycle<max_tries) and  ((frac_completeH<max_frac_completeH) or (frac_completeF<1)) ):

        print()
        #print("\rCycle: %s\t CompleteH: %s-%s,\tCompleteF: %s-%s\t%s-%s"%(cycle,frac_completeH,[x for x in H_bins_countH],frac_completeF,[x for x in H_bins_countF],sum(H_bins_countH),sum(H_bins_countF)),end="")
        for i in range(int(numEs*nbins)+1):
            # if((sum(H_bins_completedF)-1>last_sum) and (i>numEs)):
            #     last_sum=sum(H_bins_completedF)
            #     break
            print("\rCycle: %s i:%s\t%s\tH: %s-%s,\tF: %s-%s\t%s-%s\t%s\t%s"
            %(cycle,i,tshannon,frac_completeH, H_bins_countH,frac_completeF,H_bins_countF,
              sum(H_bins_countH),sum(H_bins_countF),fA,rho),end="")    
            
            this_ev=getRandomEnvironment(ndays, p=rho, alphabet=[0,1])
            
            this_ev_inv=this_ev.copy()
            this_ev_inv['Ev']=[x for x in invert_Ev(np.array(this_ev['Ev']))]
            evs=[this_ev_inv,this_ev]
            
            for this_E in evs:
                thisEv=this_E['Ev']
                fA=round(sum(thisEv)/ndays,3)
                this_E['frac_1']=fA

                tshannon=round(this_E['shannon'],4)
                ibinH=np.argmax(H_bins[:]>=tshannon)
                this_bin_countH=H_bins_countH[ibinH]

                ibinF=np.argmax(H_bins[:]>=fA)
                this_bin_countF=H_bins_countF[ibinF]
                
                if((this_bin_countH<numEs)&(this_bin_countF<numEs)):
                    Es.append(this_E)
                    Es_byBins0H[ibinH].append(env_index)
                    Es_byBins0F[ibinF].append(env_index)
                    H_bins_countH[ibinH]=this_bin_countH+1
                    H_bins_countF[ibinF]=this_bin_countF+1
                    env_index+=1
                elif((this_bin_countH<numEs)&(this_bin_countF>=numEs)):
                    Es.append(this_E)
                    Es_byBins0H[ibinH].append(env_index)
                    H_bins_countH[ibinH]=this_bin_countH+1
                    env_index+=1
                    if((this_bin_countF>=numEs)):
                        H_bins_completedF[ibinF]=1
                        frac_completeF=round(sum(H_bins_completedF)/len(H_bins_completedF),2)
                elif((this_bin_countH>=numEs)&(this_bin_countF<numEs)):
                    Es.append(this_E)
                    Es_byBins0F[ibinF].append(env_index)
                    H_bins_countF[ibinF]=this_bin_countF+1
                    env_index+=1
                    if((this_bin_countH>=numEs)):
                        H_bins_completedH[ibinH]=1
                        frac_completeH=round(sum(H_bins_completedH)/len(H_bins_completedH),2)
                elif((this_bin_countH>=numEs)&(this_bin_countF>=numEs)):
                    H_bins_completedH[ibinH]=1
                    frac_completeH=round(sum(H_bins_completedH)/len(H_bins_completedH),2)
                    H_bins_completedF[ibinF]=1
                    frac_completeF=round(sum(H_bins_completedF)/len(H_bins_completedF),2)

        cycle+=1
        rho=rho*.9
        if(rho<1e-5):
            rho=.75

    print()        
    return Es,Es_byBins0H,Es_byBins0F


def purge_repeated_envs(Es,numBins,numEs,Es_BinsH,Es_BinsF):
    newEs=[]
    H_bins=np.arange(0,1+1/numBins,1/numBins)
    H_bins_countF=[int(x) for x in np.zeros(len(H_bins))]
    all_used=[]
    for ib,this_Es_byBin in enumerate(Es_BinsH):
        for ie,this_Eindex in enumerate(this_Es_byBin):
            all_used.append(this_Eindex)
            thisE=Es[this_Eindex].copy()
            newEs.append(thisE)

            fA=thisE["frac_1"]
            ibinF=np.argmax(H_bins[:]>=fA)
            this_bin_count=H_bins_countF[ibinF]
            if this_bin_count<numEs:
                H_bins_countF[ibinF]=this_bin_count+1
    
    for ib,this_count in enumerate(H_bins_countF):
        resting=numEs-this_count
        #print(len(Es_BinsF[ib]))  
        pre_list=Es_BinsF[ib]
        choise_list=[x for x in pre_list if x not in all_used]

        if(len(choise_list)<resting):
            print("Kakatua")
            
        new_chosen=random.sample(choise_list,resting)
        #print(len(new_chosen))
        for new_i in new_chosen:   
            thisE=Es[new_i].copy()
            newEs.append(thisE)
            
    return newEs

