# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:41:13 2024

@author: mkost
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from Parameters_LV import *


'''
Returns heatmap plots. Therefore, it first calculates the particles density and then 
averages over all simulations.
'''

h=0.1

discret=int(a/h)

ts=0.1

dx_hist=a/l_coupling


sim=200

dx_hist=a/l_coupling

print(deltat)
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
    PreySimulation=list of all simulations, see Coupling.py
    '''

    
    timesteps_cut = int(deltat * (timesteps-1) / ts)
    
    number_bins=int(a/h)
    
    DiscretePrey=np.empty([sim,timesteps_cut, number_bins,number_bins])
    DiscretePreyAverage=np.empty([timesteps_cut, number_bins,number_bins])

    DiscretePred=np.empty([sim,timesteps_cut, number_bins,number_bins])
    DiscretePredAverage=np.empty([timesteps_cut, number_bins,number_bins])
    

    for s in range(sim):
       
        print(s)
        for t in range(timesteps_cut):
            
            
            Prey_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/'+str('LVPreyParticles')+str(s)+ 'time'
                    +str(t)+'.npy', allow_pickle=True)
            Predator_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/'+str('LVPredatorParticles')+str(s)+ 'time'
                    +str(t)+'.npy', allow_pickle=True)
           
            discrete_Prey_s_t=Discretization(a, l_coupling, Prey_s_t)
            
            DiscretePrey[s,t,:,:]=discrete_Prey_s_t


            discrete_Predator_s_t=Discretization(a, l_coupling, Predator_s_t)

            DiscretePred[s,t,:,:]=discrete_Predator_s_t

            
    # average over simulations 
    for t in range(timesteps_cut):
        
        DiscretePreyAverage[t,:,:]=np.mean(DiscretePrey[:,t,:,:], axis=0)
        DiscretePredAverage[t,:,:]=np.mean(DiscretePred[:,t,:,:], axis=0)   
  
    return DiscretePreyAverage, DiscretePredAverage


DiscretePreyAverage, DiscretePredAverage=functionAverage()

np.save('./Solutions/LVPreyAverage.npy', DiscretePreyAverage)
np.save('./Solutions/LVPredAverage.npy', DiscretePredAverage)


print('done averaging')
