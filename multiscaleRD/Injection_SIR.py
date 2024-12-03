# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:12:19 2020

@author: Margarita

The following functions are made to inject particles from a conctinuous domain into the particle-based domain with rate gamma.
It requires the boundaryconcentration of the boundary of width deltar and length L (the domain length).
The functions will be used in the Strang splitting in the main file (Reaction+Injection, Diffusion, Reaction+Injection)
To guarantee an accurate result, we let the particles in the boundary alo proliferate (virual proliferation).
However this particles are not allowed to be injected in the same (Strang splitting) time-step (concentrationmovement0). On the other hand the particles
that had proliferated from the PREVIOUS time-step can be injected (concentrationmovement1).

"""
import numpy as np



def concentrationmovement( Boundaryconcentration_t, deltat,deltar, L, gamma): 
    
    '''
    Returns a list of positions (2D arrays) of the new particles in the PBS domain (Children). The particles are injected with proabibility gamma. 
    Only the particles that did not proliferated in the same time-step can be injected. Therefore we have to subtract the 'Extra' particles 
    from the total boundaryconcentrations.
    deltat=time
    Booundaryconcentration_t=list or array of Boundaryconcentration of each cell of length and width deltar (L/deltar= number of boundary cells).
    deltat=time-step size
    gamme=injection rate
    deltar=boundary cell length and width
    Extra=number of  actual dead virtal preys in the boundary cell
    L=domain size
    '''
    Children=[]
    Pr=1-np.exp(-gamma*deltat) # probability of injection
    for i in range(len(Boundaryconcentration_t)):       
        
        integ, dec = int(np.floor(Boundaryconcentration_t[i])), Boundaryconcentration_t[i]-int(np.floor(Boundaryconcentration_t[i]))
        for v in range(integ): # test for every particle           
            if Pr > np.random.rand():
                Children.append(np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))]))        
        if 1-np.exp(-gamma*deltat*dec)>np.random.rand(): # for the virtual particle
            Children.append(np.array([np.random.uniform(L-deltar,L), np.random.uniform(deltar*i, deltar*(i+1))]))

    return Children

