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
import time
import multiprocessing
from multiprocessing import Pool
from functools import partial
from Reaction_SIRImplicit import *
from Injection_SIRTau import *
from ParametersB_SIR import * 



import sys

if len(sys.argv) < 2:
    print('No input provided...')
    sys.exit()
else:
    start = int(sys.argv[1]) - 1
    print('\nHello! This is task number', start)
np.random.seed(start)


'''
This main file couples a PBS to a mean-field concentration and returns trajectories of particles for different simulations. 
Before running this file it is necessary to calculate the continuous solution with FD_*.py. The functions for reactions are from Reaction_*.py and 
the one for injection from Injection.py
The input is:
D1, D2, D3 = Diffusion coefficients of S, I, R
timesteps = number of time steps
deltat = time-step size
l = number of cells in the FD scheme
a = vertical length of the domain
L = horizontal length of the domain
r1_macro = infection rate
r2 = second order microscopic rate A+B -> 2B
sigma = reaction radius
r3 = recovery rate
deltar = width of the boundary domain (one boundary cell has the length and width of deltar)
The code consists of the following components:
1) Calculate the boundary concentration for EACH timestep
2) Iteration with Strang Splitting using the function from Reaction.py and Injection.py: Injection, Reaction, Diffusion, Reaction, Injection.
3) Multiprocessing that does many simulations at the same time
'''
#%%

# Calculate the boundary concentration 'Boundaryconcentration' from the FD solution

x0 = np.array([Lx + 1, Lx])  # location of source
dx_hist = a / l_coupling  # length of histogram cell edge (squares of dx_hist by dx_hist)
yarray1 = np.arange(0, a, deltar1)  # Array to locate boundary cells for species S
yarray2 = np.arange(0, a, deltar2)  # Array to locate boundary cells for species I
yarray3 = np.arange(0, a, deltar3)  # Array to locate boundary cells for species R

# Simulation parameters
dtshould1 = deltar1 * deltar1 / (2.0 * D1)  # time-step size
dtshould2 = deltar2 * deltar2 / (2.0 * D2)  # time-step size
dtshould3 = deltar3 * deltar3 / (2.0 * D3)  # time-step size
print(dtshould1 , 'dtmax', deltar1,'dr1', deltar2,'dr2', deltar3, 'dr3', r1,'r1', r2,'r2', sigma, 'sigma')


# load concentration

listC1 = np.load('./Solutions/FDSIR1_B.npy')  # gets Data from continuous solution of S
listC2 = np.load('./Solutions/FDSIR2_B.npy')  # gets Data from continuous solution of I
listC3 = np.load('./Solutions/FDSIR3_B.npy')  # gets Data from continuous solution of R


# transalte into number of particles in each boundary cell 

averageNumberParticles1 = np.zeros((len(yarray1), timesteps))
averageNumberParticles2 = np.zeros((len(yarray2), timesteps))
averageNumberParticles3 = np.zeros((len(yarray3), timesteps))

xlimits1 = [a / 2, a / 2 + deltar1]
xlimits2 = [a / 2, a / 2 + deltar2]
xlimits3 = [a / 2, a / 2 + deltar3]

gamma1 = D1 / ((deltar1) ** 2)  # injection rate
gamma2 = D2 / ((deltar2) ** 2)
gamma3 = D3 / ((deltar3) ** 2)

# Function to calculate number of particles in each boundary
def calculate_boundary_concentration(yarray, deltar, listC, gamma):
    averageNumberParticles = np.zeros((len(yarray), timesteps))
    for i in range(len(yarray)):
        ylimits = [yarray[i], yarray[i] + deltar]
        for k in range(timesteps):
            if k == 0:
                averageNumberParticles[i, k] = 0.0
            else:
                averageNumberParticles[i, k] = deltar * deltar * listC[k][int(i / (len(yarray) / l_coupling)), int(l_coupling / 2)]
    return averageNumberParticles

Boundaryconcentration1 = calculate_boundary_concentration(yarray1, deltar1, listC1, gamma1)
Boundaryconcentration2 = calculate_boundary_concentration(yarray2, deltar2, listC2, gamma2)
Boundaryconcentration3 = calculate_boundary_concentration(yarray3, deltar3, listC3, gamma3)
# Iteration
#%%
start_time=time.time()

def functionsimulation(ts):
    internal_reaction_total_time=0
    '''
    Returns lists for S, I, R consisting of trajectories at every desirable time-step and for every simulation and the reference solutions
    simulations = number of simulations
    ts = which timestep we save: 0, ts, ts*2
    Sus = S
    Inf = I
    Recovery = R
    Children refers to newly infected people (the name comes from LV dynamics)
    SusSimulation = saves all simulations
    SusPosition = updating list containing the CURRENT position of susceptible people as a 
    2D array
    SusPositionHalfTime = contains the positions at each desirable time step
    Structure: SusSimulation = [[Simulation_1], [Simulation_2]...], Simulation_1 = [[Time_1], [Time_2],..], for example Time_1 = [[1,2], [3,4],[1.4,4]] (positions)
    Analogous for infected
    '''

    SusPosition = []
    RecoveryPosition = []
    InfPosition = []
    
    timesteps_cut = int(deltat * (timesteps-1) / ts)
    
    SusPositionHalfTime = [[] for _ in range(timesteps_cut)]
    InfPositionHalfTime = [[] for _ in range(timesteps_cut)]
    RecoveryPositionHalfTime = [[] for _ in range(timesteps_cut)]
    
    k = 0
    
    for t in range(timesteps):
        
        # Injection
       
        SusChildrenInj = concentrationmovement_tauleaping(Boundaryconcentration1[:, t], deltat * 1 / 2, deltar1, Lx,
                                                  gamma1,1)
        InfChildrenInj = concentrationmovement_tauleaping(Boundaryconcentration2[:, t], deltat * 1 / 2, deltar2, Lx,
                                                 gamma2,1)
        RecoveryChildrenInj = concentrationmovement_tauleaping(Boundaryconcentration3[:, t],deltat * 1 / 2, deltar3, Lx
                                                    , gamma3,1)

        # Reaction
        # Recovery of I
        InfPosition, RecoveryNew1 = dying(InfPosition, r2, deltat * 0.25, InfPosition)
        # Reactions between S and I
        SusPosition, InfChildrenReact, InfB, internal_time = eatcompact(SusPosition, InfPosition, Lx, deltar1, deltar2,
                                                            Boundaryconcentration1[:, t],
                                                            Boundaryconcentration2[:, t], r1, sigma, deltat*0.5,1)
        
        internal_reaction_total_time += internal_time

        # Recovery of Infected
        InfPosition, RecoveryNew2 = dying(InfPosition, r2, deltat * 0.25, InfB)

        # Put them all together, such that new injected particles can react
        SusPosition.extend(SusChildrenInj)
        InfPosition.extend(InfChildrenInj + InfChildrenReact)
        RecoveryPosition.extend(RecoveryChildrenInj + RecoveryNew1 + RecoveryNew2)

        # Diffusion
        SusPosition = movement(SusPosition, deltat, D1, Lx, a)
        InfPosition = movement(InfPosition, deltat, D2, Lx, a)
        RecoveryPosition = movement(RecoveryPosition, deltat, D3, Lx, a)

        # Reaction
        InfPosition, RecoveryNew1 = dying(InfPosition, r2, deltat * 0.25, InfPosition)
        SusPosition, InfChildrenReact, InfB, internal_time = eatcompact(SusPosition, InfPosition, Lx, deltar1, deltar2,
                                                            Boundaryconcentration1[:, t],
                                                            Boundaryconcentration2[:, t], r1, sigma, deltat*0.5, 1)
        
        internal_reaction_total_time += internal_time
        InfPosition, RecoveryNew2 = dying(InfPosition, r2, deltat * 0.25, InfB)

        # Injection
        SusChildrenInj = concentrationmovement_tauleaping(Boundaryconcentration1[:, t], deltat * 1 / 2, deltar1, Lx,
                                                  gamma1,1)
        InfChildrenInj = concentrationmovement_tauleaping(Boundaryconcentration2[:, t], deltat * 1 / 2, deltar2, Lx,
                                                 gamma2,1)
        RecoveryChildrenInj = concentrationmovement_tauleaping(Boundaryconcentration3[:, t], deltat * 1 / 2, deltar3, Lx
                                                    , gamma3,1)

        # Put them all together
        SusPosition.extend(SusChildrenInj)
        InfPosition.extend(InfChildrenInj + InfChildrenReact)
        RecoveryPosition.extend(RecoveryChildrenInj + RecoveryNew1 + RecoveryNew2)

        if t%int(timesteps/timesteps_cut)==0 and t!=0:
            
            SusPositionHalfTime[k] = np.array(SusPosition)
            InfPositionHalfTime[k] = np.array(InfPosition)
            RecoveryPositionHalfTime[k] = np.array(RecoveryPosition)
            # save results
            np.save('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIRSusTauBNew' + str(start) + 'time'
                    +str(k)+'', SusPositionHalfTime[k])
            np.save('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIRInfTauBNew' + str(start) + 'time'
                    +str(k)+'',
                    InfPositionHalfTime[k])
            np.save('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIRRecoveryTauBNew' + str(start) +'time'
                    +str(k)+ '',
                    RecoveryPositionHalfTime[k])
            k += 1
            print(t)
    
    
    return SusPositionHalfTime, InfPositionHalfTime, RecoveryPositionHalfTime,internal_reaction_total_time
#%%
save_time=0.1

# Call the function
SusSimulation, InfSimulation,  RecoverySimulation,internal_reaction_total_time= functionsimulation(save_time)
end_time=time.time()

print('Task number', start, 'is done')
# Calculate elapsed time
elapsed_time = end_time - start_time-internal_reaction_total_time

# Save the elapsed time as a .npy file
np.save('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIR_timeTauCutBNew'+str(start)+'.npy', elapsed_time)

