#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 14:20:10 2020

@author: bzfkostr
"""
from __future__ import division
import numpy as np

'Parameters'

def fconstant(x):
	if 5.1<=x[0]<=7.1 and 4<=x[1]<=6 :
		ux=100
	else:
		ux=0
	
	return ux    

def fconstant3(x):
	if 5.6<=x[0]<=6.6 and 4.5<=x[1]<=5.5 :
		ux=10
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

Computati0nal Parameters
sim_number=number of simulations in each parallel simulation
numSimulations=how many simulations should run in parallel
s=which timestep we want to save

'''

uPrey=fconstant
uPred=fconstant3
l=100+1
m=2*l-1
m=l
deltat=0.0025
timesteps=1000
a=10
h=a/(l-1) # grid size 
D1=0.3
D2=0.1
r1=0.15 # micro = macro parameters
r3=0.1  # micro = macro parameters   


''' Parameters cor coupling'''

sigma=0.01
r2_macro=0.05
r2=r2_macro/(np.pi*(sigma**2)) 
L=a/2 # half of the grid
deltar1=np.sqrt(deltat*D1*2) # boundary size
deltar2=np.sqrt(deltat*D2*2) # boundary size
l_coupling=l-1 
sim_number=20 
numSimulations = 30
s=0.1 

