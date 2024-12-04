#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:23:43 2019

@author: bzfkostr
"""
from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("agg")
from Reaction_LV import *
from Injection_LV import *
from Parameters_LV import *
import sys

#code for LV parameters (not published)

# multiprocessing

if len(sys.argv) < 2:
    print('No input provided...')
    sys.exit()
else:
    start = int(sys.argv[1]) - 1
    print('\nHello! This is task number', start)

'''
This main file couples a PBS to a mean-field concentration and returns trajectories of particles for different simulations. 
Before running this file it is necessary to calculate the continuous solution with FD.py. The functions for reactions are from Reaction.py and 
the one for injection from Injection.py
The input is:
D1, D2=Diffusion coefficient of A and B
timesteps=stepsnumber of time
deltat=time-step size
l=number of cells in the FD scheme
a=vertical length of the domain
L=horizontal length of the domain
r1=first order microscopic rate A-> 2A
r2_macro=second order macroscopic rate, that corresponds to the rate in the FD scheme
r2=second order microscopic rate A+B-> 2B
sigma=reaction radius
r3=zero reaction B->0 microscopic rate
deltar=width of the boundary domain (one boundary cell has the length and width of deltar)
The code consists of the following components
1) Calculate the boundary concentration for EACH timestep
2) Iteration with Strang Splitting using the function from Reaction.py and Injection.py: Injection, Reaction, Diffusion, Reaction, Injection.
3) Multiprocessing that does many simulations at same time
'''

# Calculate the boundary concentration 'Boundaryconcentration' from the FD solution

x0 = np.array([L+1, L])  # location of source
dx_hist = a / l_coupling  # length of histogram cell edge (squares of dx_hist by dx_hist)
yarray1 = np.arange(0, a, deltar1)  # Array to locate boundary cells for species A
yarray2 = np.arange(0, a, deltar2)  # Array to locate boundary cells for species B

# Simulation parameters
dtshould1 = deltar1 * deltar1 / (2.0 * D1)  # time-step size
dtshould2 = deltar2 * deltar2 / (2.0 * D2)  # time-step size
print(dtshould1, dtshould2, deltat, 'Should be equal')

#Time = np.linspace(0, maxtime, timesteps)

listC1 = np.load('./Solutions/FDLVSolution1.npy')  # gets Data from continuous solution of A
listC2 = np.load('./Solutions/FDLVSolution2.npy')  # gets Data from continuous solution of B

averageNumberParticles1 = np.zeros((len(yarray1), timesteps))
averageNumberParticles2 = np.zeros((len(yarray2), timesteps))
xlimits1 = [a / 2, a / 2 + deltar1]
xlimits2 = [a / 2, a / 2 + deltar2]

gamma1 = D1 / ((deltar1) ** 2)  # injection rate
gamma2 = D2 / ((deltar2) ** 2)

# Boundary concentration for A

for i in range(len(yarray1)):
    ylimits1 = [yarray1[i], yarray1[i] + deltar1]

    for k in range(timesteps):
        if k == 0:
            averageNumberParticles1[i, k] = 0.0
        else:
            averageNumberParticles1[i, k] = (deltar1) * (deltar1) * listC1[k][int(i / (len(yarray1) / l_coupling)), int(
                l_coupling / 2)]

# Boundary concentration for B

for i in range(len(yarray2)):
    ylimits2 = [yarray2[i], yarray2[i] + deltar2]
    for k in range(timesteps):
        if k == 0:
            averageNumberParticles2[i, k] = 0.0
        else:
            averageNumberParticles2[i, k] = (deltar2) * (deltar2) * listC2[k][int(i / (len(yarray2) / l_coupling)),
                                                                           int(l_coupling / 2)]

Boundaryconcentration1 = averageNumberParticles1
Boundaryconcentration2 = averageNumberParticles2

# Iteration
def arrays_to_tuples(arrays):
    return set(tuple(array) for array in arrays)

def functionsimulation(ts):
    # which time-steps to save
    timesteps_cut = int(deltat * (timesteps-1) / ts)
    
    np.random.seed(int(start))

    PreyPosition = []
    PredatorPosition = []
    
    PreyPositionHalfTime = [[] for _ in range(timesteps_cut)]
    PredatorPositionHalfTime = [[] for _ in range(timesteps_cut)]

    k = 0
    for t in range(timesteps):

        # Injection
        PreyChildrenInj = concentrationmovement(Boundaryconcentration1[:, t], deltat * 1 / 2, deltar1, L, gamma1)
        PredChildrenInj = concentrationmovement(Boundaryconcentration2[:, t], deltat * 1 / 2, deltar2, L, gamma2)

        # Reaction
        PreyChildrenProlif1, NotProliferatedPreys = proliferation(PreyPosition, r1, deltat * 0.25)
        PredatorPosition = dying(PredatorPosition, r3, deltat * 0.25, PredatorPosition)
        PreyPositionAfterReaction, PredChildrenReact, PredB = eatcompact(NotProliferatedPreys, PredatorPosition, L, deltar1,deltar2,
                                                            Boundaryconcentration1[:, t],
                                                            Boundaryconcentration2[:, t], r2, sigma, deltat*0.5, L)
       
        set_A = arrays_to_tuples(PreyPositionAfterReaction)
        set_B =  arrays_to_tuples(NotProliferatedPreys)
        difference = set_B - set_A
        NotReactedPreys = [np.array(list(tuples)) for tuples in difference]

        PreyChildrenProlif2, NotProliferatedPreys = proliferation(NotReactedPreys, r1, deltat * 0.25)
        PredatorPosition = dying(PredatorPosition, r3, deltat * 0.25, PredB)

        # Put them all together
        PreyPosition = PreyChildrenInj + PreyChildrenProlif1 + PreyPositionAfterReaction + PreyChildrenProlif2
        PredatorPosition = PredChildrenInj + PredatorPosition + PredChildrenReact

        # Diffusion
        PreyPosition = movement(PreyPosition, deltat, D1, L, a)
        PredatorPosition = movement(PredatorPosition, deltat, D2, L, a)

        # Reaction
        PreyChildrenProlif1, NotProliferatedPreys = proliferation(PreyPosition, r1, deltat * 0.25)
        PredatorPosition = dying(PredatorPosition, r3, deltat * 0.25, PredatorPosition)
        PreyPositionAfterReaction, PredChildrenReact, PredB = eatcompact(NotProliferatedPreys, PredatorPosition, L, deltar1,deltar2,
                                                            Boundaryconcentration1[:, t],
                                                            Boundaryconcentration2[:, t], r2, sigma, deltat*0.5, L)
        
        
        set_A = arrays_to_tuples(PreyPositionAfterReaction)
        set_B = arrays_to_tuples(NotProliferatedPreys)
        difference = set_B - set_A
        NotReactedPreys = [np.array(list(tuples)) for tuples in difference]

        PreyChildrenProlif2, NotProliferatedPreys = proliferation(NotReactedPreys, r1, deltat * 0.25)
        PredatorPosition = dying(PredatorPosition, r3, deltat * 0.25, PredB)

        PreyChildrenInj = concentrationmovement(Boundaryconcentration1[:, t], deltat * 1 / 2, deltar1, L, gamma1)
        PredChildrenInj = concentrationmovement(Boundaryconcentration2[:, t], deltat * 1 / 2, deltar2, L, gamma2)

        # Put them all together
        PreyPosition = PreyChildrenInj + PreyChildrenProlif1 + PreyPositionAfterReaction + PreyChildrenProlif2
        PredatorPosition = PredChildrenInj + PredatorPosition + PredChildrenReact

    
       
        if t%int((timesteps-1)/timesteps_cut)==0 and t!=0:
                
            PreyPositionHalfTime[k] = np.array(PreyPosition)
            PredatorPositionHalfTime[k] = np.array(PredatorPosition)
            
            np.save(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/LVPreyParticles{start}time{k}', PreyPositionHalfTime[k])
            np.save(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/LVPredatorParticles{start}time{k}', PredatorPositionHalfTime[k])
            
            k += 1
            print((t)*deltat, timesteps/timesteps_cut)
    
    return PreyPositionHalfTime, PredatorPositionHalfTime
#%%
save_time=0.1

# Call the function
PreySimulation, PredatorSimulation= functionsimulation(save_time)
