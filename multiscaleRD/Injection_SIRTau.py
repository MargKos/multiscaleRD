# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:12:19 2020

@author: Margarita

The following functions are made to inject particles from a conctinuous domain into the particle-based domain with rate gamma.
It requires the boundaryconcentration of the boundary of width deltar and length L (the domain length).
The functions will be used in the Strang splitting in the main file (Reaction+Injection, Diffusion, Reaction+Injection)
"""


import numpy as np
def concentrationmovement_tauleaping(Boundaryconcentration_t, tau, deltar, L, gamma, M): 
    '''
    Returns a list of positions (2D arrays) of the new particles in the PBS domain (Children). 
    The particles are injected with probability gamma. 

    Parameters:
    - Boundaryconcentration_t: list or array of Boundary concentration of each cell of length and width deltar (L/deltar= number of boundary cells).
    - M: number of sub-steps in tau-leaping.
    - tau: overall time step size.
    - deltar: spatial cell size.
    - L: domain size.
    - gamma: injection rate.
    
    Returns:
    - Children: List of 2D arrays with positions of newly injected particles in the PBS domain.
    '''
    
    Children = []
    
    # Loop over all boundary cells
    for i in range(len(Boundaryconcentration_t)):    
        
        particles_in_cell = Boundaryconcentration_t[i]  # Number of particles in cell i
        injected_particles = 0  # Number of injected particles initially 0
        delta_tau = tau / M  # Time step for each sub-step in tau-leaping
    
        # Loop over M sub-steps
        Lambda = gamma * (particles_in_cell - injected_particles)
        for k in range(M):
            if particles_in_cell <= injected_particles:
                break  # Exit loop if no more particles to inject
            
            
            injected_particles += np.random.poisson(Lambda * delta_tau)  # Update injected_particles

        # Inject particles into the domain based on the number injected
        for j in range(injected_particles):
            x_pos = np.random.uniform(L - deltar, L)
            y_pos = np.random.uniform(deltar * i, deltar * (i + 1))
            Children.append(np.array([x_pos, y_pos]))
        
        

    return Children
