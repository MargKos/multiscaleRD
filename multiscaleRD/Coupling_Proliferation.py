#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:23:43 2019

@author: bzfkostr
"""

from __future__ import division
import numpy as np
import matplotlib
import time
matplotlib.use("agg")
from functools import partial
from Reaction_LV import *
from Injection_LV import *
from Parameters_Proliferation import *
import sys


if len(sys.argv) < 2:
    print('No input provided...')
    sys.exit()
else:
    start = int(sys.argv[1]) - 1
    print('\nHello! This is task number', start)

'''
This main file couples a PBS to a mean-field concentration and returns trajectories of particles for different simulations. 
Before running this file, it is necessary to calculate the continuous solution with FD_Proliferation.py. The functions for reactions are from Reaction_LV.py and 
the one for injection from Injection_LV.py.

Input parameters:
- D: Diffusion coefficient
- timesteps: Number of time steps
- deltat: Time-step size
- l: Number of cells in the FD scheme
- a: Vertical length of the domain
- L: Horizontal length of the domain
- r1: First-order microscopic rate
- deltar: Width of the boundary domain (one boundary cell has the length and width of deltar)

The code consists of the following components:
1) Calculate the boundary concentration for EACH timestep
2) Iteration with Strang Splitting using functions from Reaction.py and Injection.py: Injection, Reaction, Diffusion, Reaction, Injection
'''


# multiprocessing


deltar = np.sqrt(deltat * D * 2) 

x0 = np.array([L + 1, L])  # Location of source
dx_hist = a / l_coupling  # Length of histogram cell edge (equals to h)
dtshould = deltar * deltar / (2.0 * D)  # Time-step size
print(dtshould, deltat, 'Should be equal')
gamma = D / (deltar ** 2)  # Injection rate

'''1. Calculate the boundary concentration 'Boundaryconcentration' from the FD solution'''

maxtime = deltat * (timesteps - 1)  # Maximum time simulation can reach
Time = np.linspace(0, maxtime, timesteps)
listC = np.load('./Solutions/FDSolution_Proliferation.npy')  # Data from continuous solution
yarray = np.arange(0, a, deltar)  # Array to locate boundary cells  
averageNumberParticles = np.zeros((len(yarray), timesteps))

for i in range(len(yarray)):
    ylimits = [yarray[i], yarray[i] + deltar]
    for k in range(timesteps):
        if k == 0:
            averageNumberParticles[i, k] = 0.0
        else:
            averageNumberParticles[i, k] = (deltar ** 2) * listC[k][int(i / (len(yarray) / l_coupling)), int(l_coupling / 2)]  # X = C * V

Boundaryconcentration = averageNumberParticles

# Iteration
def arrays_to_tuples(arrays):
    return set(tuple(array) for array in arrays)
    
start_time = time.time()
def functionsimulation(ts):
    # which time-steps to save
    timesteps_cut = int(deltat * (timesteps-1) / ts)
    
    np.random.seed(int(start))

    PreyPosition = []
   
    
    PreyPositionHalfTime = [[] for _ in range(timesteps_cut)]
   
    print(np.shape(PreyPositionHalfTime), ts)
    k = 0
    for t in range(timesteps):

        # Injection
        PreyChildrenInj = concentrationmovement(Boundaryconcentration[:, t], deltat * 0.5, deltar, L, gamma)
    

        # Reaction
        PreyChildrenProlif1, NotProliferatedPreys = proliferation(PreyPosition, r1, deltat * 0.5)
       

        # Put them all together
        PreyPosition = PreyChildrenInj + PreyChildrenProlif1 + PreyPosition


        # Diffusion
        PreyPosition = movement(PreyPosition, deltat, D, L, a)
      
        # Reaction
        PreyChildrenProlif1, NotProliferatedPreys = proliferation(PreyPosition, r1, deltat * 0.5)
       
        PreyChildrenInj = concentrationmovement(Boundaryconcentration[:, t], deltat * 0.5, deltar, L, gamma)


        # Put them all together
        PreyPosition = PreyChildrenInj + PreyChildrenProlif1 + PreyPosition
       
        if t%int((timesteps-1)/timesteps_cut)==0 and t!=0:
                
            PreyPositionHalfTime[k] = np.array(PreyPosition)
            np.save(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/ProliferationParticles_{start}_time{k}.npy', PreyPosition)
            
            k += 1
            print((t)*deltat, timesteps/timesteps_cut)
    
    return PreyPositionHalfTime
#%%
save_time=0.1

# Call the function
PreySimulation= functionsimulation(save_time)
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time

# Save the elapsed time as a .npy file
np.save('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/execution_time'+str(start)+'.npy', elapsed_time)