#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 15:46:15 2020

@author: bzfkostr
"""

import numpy as np

c=250000
def u_0(x):
	if 5<x[0]<=5.5 :
		ux=c*(2*abs(np.sin(2*np.pi*x[1]/10)))
	else:
		ux=0
	
	return ux    
'Parameters'

'''
Space discretization
l+1=number of grid-cells in y direction
m+1=number of grid-cells in x direction

Time discretization
deltat=time-step size
timesteps=number of iterations 

Mathematical parameters
D=diffusion coefficient
r=first order macroscopic rate
a=domain length
l_coupling=number of cells in the FD scheme after the averaging
a=vetical length of the doman
L=horizonal length of the domain equals to the half of a
deltar=width of the boundary domain (one boundary cell has the legnth and width of deltar)

Computational Parameters
sim_number=number of simulations in each parallel simulation
numSimulations=how many simulations should run in parallel
s=which timestep we want to save
'''


l=100+1 # number of grids in x direction, should fit to a s.t. a/(l-1) is a 'good' number
m=l # number of grids in y direction

deltat=0.002
timesteps=500+1 

D=3
a=10 #domain boundary

  

''' Parameters only for Coupling'''
L=a/2 # x coordinate where the the PBS and FD domain are seperated
deltar=np.sqrt(deltat*D*2) # boundary size
l_coupling=l-1 
h=0.1
#sim_number=500

