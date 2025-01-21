# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:41:13 2024

@author: mkost
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from ParametersB_SIR import *


'''
Returns heatmap plots. Therefore, it first calculates the particles density and then 
averages over all simulations.
'''

h=0.1

discret=int(a/h)

ts=0.1

dx_hist=a/l_coupling

maxtime = deltat * (timesteps - 1)  # maximum time simulation can reach


sim=500

dx_hist=a/l_coupling

'''Densities'''

def Discretization(a, discret, Particles): 
    
    '''
    Return the concentration of particles.
    a=domain length
    discret=discrtization parameter (number of cells)
    Particles=list of 2D arrays
    '''
    
    xPositions=[]
    yPositions=[]
    for i in range(len(Particles)):
       
        xPositions.append(Particles[i][0])
        yPositions.append(Particles[i][1])
    
   
    xbins=np.arange(0,a+h, h)
    ybins=np.arange(0,a+h, h)
    concentration, xbins, ybins=np.histogram2d(xPositions, yPositions, bins=(xbins, ybins), weights=np.ones_like(xPositions)/(h**2))
    
    return concentration
    

def functionAverage(): 
    
    '''
    Returns the mean-field concentration for each time-step by averaging over all
    simulations
    '''

    
    k=0
    
    timesteps_cut=int(deltat*(timesteps-1)/ts)
    number_bins=int(a/h)
    
    DiscreteSus=np.empty([sim,timesteps_cut, number_bins,number_bins])
    DiscreteSusAverage=np.empty([timesteps_cut, number_bins,number_bins])
    
    DiscreteRecovery=np.empty([sim, timesteps_cut, number_bins,number_bins])
    DiscreteRecoveryAverage=np.empty([timesteps_cut, number_bins,number_bins])
    
    DiscreteInf=np.empty([sim,timesteps_cut, number_bins,number_bins])
    DiscreteInfAverage=np.empty([timesteps_cut, number_bins,number_bins])
 
    for s in range(sim):

        print(s)
        for t in range(timesteps_cut):

            Sus_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscaleFinal/'+str('SIRSusTauB')+str(s)+ 'time'
                    +str(t)+'.npy', allow_pickle=True)
            Inf_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscaleFinal/'+str('SIRInfTauB')+str(s)+ 'time'
                    +str(t)+'.npy', allow_pickle=True)
            Recovery_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscaleFinal/'+str('SIRRecoveryTauB')+str(s)+ 'time'
                    +str(t)+'.npy', allow_pickle=True)
            
            # Discretize Sus
            
            discrete_Sus_s_t=Discretization(a, l_coupling, Sus_s_t)
            
            DiscreteSus[s,t,:,:]=discrete_Sus_s_t

            # Discretize Inf

            discrete_Inf_s_t=Discretization(a, l_coupling, Inf_s_t)

            DiscreteInf[s,t,:,:]=discrete_Inf_s_t

            # Discretize Rec


            discrete_Recovery_s_t=Discretization(a, l_coupling, Recovery_s_t)

            DiscreteRecovery[s,t,:,:]=discrete_Recovery_s_t

    
    # average over simulations 
    for t in range(timesteps_cut):
        
        DiscreteSusAverage[t,:,:]=np.mean(DiscreteSus[:,t,:,:], axis=0)
        DiscreteInfAverage[t,:,:]=np.mean(DiscreteInf[:,t,:,:], axis=0)
        DiscreteRecoveryAverage[t,:,:]=np.mean(DiscreteRecovery[:,t,:,:], axis=0)
           

    return DiscreteSusAverage, DiscreteInfAverage, DiscreteRecoveryAverage


DiscreteSusAverage, DiscreteInfAverage, DiscreteRecAverage=functionAverage()

np.save('./Solutions/SusAverageTauB.npy', DiscreteSusAverage)
np.save('./Solutions/InfAverageTauB.npy', DiscreteInfAverage)
np.save('./Solutions/RecAverageTauB.npy', DiscreteRecAverage)

print('done averaging')
