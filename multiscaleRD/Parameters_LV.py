#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 14:20:10 2020

@author: bzfkostr
"""
from __future__ import division
import numpy as np

'Parameters for LV with high concentrations and same diffusion coefficient'

def fconstant1(x):
	if 5.1<=x[0]<=5.6 and 3.5<=x[1]<=6.5 :
		ux=2000
	else:
		ux=0
	
	return ux    

def fconstant2(x):
	if 6.1<=x[0]<=6.6 and 4.5<=x[1]<=5.5 :
		ux=1300
	else:
		ux=0
	
	return ux 
'''
Space discretization
l+1=number of grid-cells in y direction
m+1=number of grid-cells in x direction

Time discretization
deltat=time-step size
n=number of iterations 

Mathematical parameters
D1, D2=diffusion coefficients of u1 and u2
r1=first order macro and microscopic rate
r2_macro=second oder macroscopic rate
r2= microscopic rate
r3=zero reaction macro and microscopic rate
L=domain length

REMINDER: Lotka-Volterra equation
A->2A with rate r1
A+B->2B with rate r2_macro
B->0 with rate r3

A=Preys
B=Predators

Computational Parameters
sim_number=number of simulations in each parallel simulation
numSimulations=how many simulations should run in parallel
s=which timestep we want to save

'''

uPrey=fconstant1 # set initial conditions for Preys
uPred=fconstant2 # set initial conditions for Predators
l=100+1
m=l
deltat=0.0025
timesteps=1000+1
maxtime=(timesteps-1)*deltat
a=10
h=a/(l-1) # grid size 
D1=0.3
D2=0.3
r1=0.15 # micro = macro parameters, proliferation
r3=0.3  # micro = macro parameters, dying   


''' Parameters for coupling'''

sigma=0.01
r2_macro=0.01
r2=r2_macro/(np.pi*(sigma**2)) 
L=a/2 # half of the grid
deltar1=np.sqrt(deltat*D1*2) # boundary size
deltar2=np.sqrt(deltat*D2*2) # boundary size
l_coupling=l-1 
save_time=0.1 