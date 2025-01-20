# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:42:00 2020

@author: Margarita
"""

from __future__ import division
import numpy as np
from math import exp
import time

# in this code prey and predator have the same size of boundary cells

def movement(Particles,  deltat, D, Lx, a):
    
    ''' 
    Returns the new position of particles as a list after one time-step of length deltat:
    Particles=list of 2D arrays
    deltat=time-step size
    D= diffusion coefficient
    Lx=length of the horozontal boundary
    Ly=length of the vetical boundary
    We use the  Euler-Maruyama scheme, and REFLECTING boundary conditions.
    '''
    
    for i in reversed(range(len(Particles))):
        Particles[i]=np.array(Particles[i])+np.random.normal(0,np.sqrt(2*deltat*D),2) 
        if Particles[i][0]>= Lx: 
            Particles.pop(i)
        else: 
            if Particles[i][0]<=0:
               Particles[i][0]=-Particles[i][0]
            
            if Particles[i][0]>=a:
               Particles[i][0]=a-(Particles[i][0]-a)
                
            if Particles[i][1]>=a:
                Particles[i][1]=a-(Particles[i][1]-a)
            if Particles[i][1]<=0:
                Particles[i][1]=-Particles[i][1]
        
    return Particles


def euclid(v1, v2):
    
    '''
    Returns the euclidean distance between arrays v1 and v2 in an
    efficient way
    v1,v2=vectors with the same length
    '''
    dist = [(a-b)**2 for a,b in zip(v1, v2)]
    dist = np.sqrt(sum(dist))
    return dist

def dying(Particles,rate, deltat, NotImmune):
    
    '''
    Simulates the recovery
    It returns a new list of arrays (position of the species) with removed
    infected species.and a list of new recpvered species
    Particles=list of 2D arrays (positions)
    rate=microscopic reaction rate
    deltat=time-step size
    NotImmune= positition of Particles that actually are able to recover
    If all particles can die, just set NotImmune=Particles
    '''
    
    pdying=1-np.exp(-rate*deltat) 
    Recovered=[]
    for r in reversed(range(len(Particles))):
        if pdying > np.random.rand() and any((Particles[r] == x).all() for x in NotImmune):
            Recovered.append(Particles[r])
            Particles.pop(r)
            
    return Particles, Recovered

def second_order_reaction(A, B, rate, deltat, sigma, Lx):
    
    '''Simulates the second order reaction S+I->2I, if S and I are closer then sigma 
    (reaction radius). It returns the new list of S particles with removed particles, a list of 
    I particles with new particles and a list of the new particles (new list of B=previous list of B +children)
    A, B=list of 2D arrays
    A=Sus/Prey
    B=Inf/Pred
    rate=microscopic reaction rate
    deltat=time-step size
    sigma=reaction radi
    '''
  
    p=1-np.exp(-rate*deltat)
    CH=[]
    for i in reversed(range(len(A))):
        for j in reversed(range(len(B))):
            if euclid(A[i], B[j])<= sigma and p>np.random.rand():
                
                child=B[j]
               
                CH.append(child)
                    
                B.pop(j)
                A.pop(i)
                break

    return A, CH, B


def virtual(Lx, deltar, N,i):
    
    
    
    Virtual=[None]*N
    for j in range(N):
        
        Virtual[j]=np.array([np.random.uniform(Lx, Lx+deltar), np.random.uniform(deltar*i,deltar*(i+1))])
    return Virtual

def eatcompact(A, B, Lx, deltar1,deltar2,BC1, BC2, rate, sigma, deltat):
    '''
    Hybrid algorithm for second order reaction S+I->2I, if S and I are closer than sigma 
    (reaction radius). It returns the new list of S particles with removed particles, a list of 
    I particles with new particles and a list of the new particles (new list of B = previous list of B + children).
    A, B = list of 2D arrays
    deltar = length of boundary domain
    BC1 = boundary concentration of A particles
    BC2 = boundary concentration of B particles
    rate = microscopic reaction rate
    sigma = reaction radius
    L = x-boundary coordinate
    deltat = time-step size
    '''
    
    kindergarten = []

    # Create all Virtual Predators and Preys
    VirtualPreys = []
    for i in range(len(BC1)): 
        integer = int(BC1[i])
        
        VirtualPreys.extend(virtual(Lx, deltar1, integer, i))

    VirtualPreds = []
    for i in range(len(BC2)):
        integer = int(BC2[i])
        VirtualPreds.extend(virtual(Lx, deltar2, integer, i))

    # Combine actual and virtual particles for reactions
    AllPreys = A + []
    AllPreds = B + []
    
    # reactions between particles in the particle domain
    SurvivedPreys, children, NotReactedPred = second_order_reaction(AllPreys, AllPreds, rate, deltat, sigma, Lx)
    
    kindergarten.extend(children)
    
    # reactions between preys in the c domain and preds in PB
    
    PreysinCD, childrenvirtual, NotReactedPred = second_order_reaction(VirtualPreys, NotReactedPred, rate, deltat, sigma, Lx)
    
    
    # reactions between preds in the c domain and preys in PB
    internal_start = time.time()
    SurvivedPreys, children, PredsinCD = second_order_reaction(SurvivedPreys, VirtualPreds, rate, deltat, sigma, Lx)
    
    kindergarten.extend(children)
    internal_end = time.time()

    internal_time = internal_end-internal_start
    
    DecPreys = np.zeros(len(BC1))
    SingleVirtualPrey = []
    for i in range(len(BC1)):
         dec = BC1[i] - int(BC1[i])
         SingleVirtualPrey.append(virtual(Lx, deltar1, 1, i))
         DecPreys[i] = dec

    DecPreds = np.zeros(len(BC2))
    SingleVirtualPred = []
    for i in range(len(BC2)):
         dec = BC2[i] - int(BC2[i])
         SingleVirtualPred.append(virtual(Lx, deltar2, 1, i))
         DecPreds[i] = dec

    for i in range(len(SingleVirtualPrey)):
        Ignore, children_in_cdomain, NotReactedPred = second_order_reaction(SingleVirtualPrey[i], NotReactedPred, rate * DecPreys[i], deltat, sigma, Lx)
    
    for i in range(len(SingleVirtualPred)):
        SurvivedPreys, childrenvirtual,VirtualPreds = second_order_reaction(SurvivedPreys, SingleVirtualPred[i], rate * DecPreds[i], deltat, sigma, Lx)
        kindergarten.extend(childrenvirtual)
       

    return SurvivedPreys, kindergarten, NotReactedPred, internal_time
