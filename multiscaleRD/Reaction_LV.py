# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:42:00 2020

@author: Margarita
"""

from __future__ import division
import numpy as np
from math import exp
''' Code for LV open system with different diffusion coefficients of prey and predator with explicit reaction across borders'''



def movement(Particles,  deltat, D, Lx, Ly):
    
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
                
            if Particles[i][1]>=Ly:
                Particles[i][1]=Ly-(Particles[i][1]-Ly)
            if Particles[i][1]<=0:
                Particles[i][1]=-Particles[i][1]
        
    return Particles




def proliferation(Particles, rate, deltat):
    
    ''' 
    Simulates the first oder reaction A->2A or the proliferation of a species. It returns a list of positions (2D array)
    of new particles (children).
    Particles=list of 2D arrays
    rate=microscopic reaction rate of the first order reaction
    deltat=time-step size
    return positions of new particles, that are equal to the position of the proliferated particle, and position of the non proliferated particles
    '''
  
    pproliferation=1-np.exp(-rate*deltat) 
    
    children=[] 
    for r in reversed(range(len(Particles))): 
        if pproliferation > np.random.rand():
           children.append(np.array(Particles[r]))
           Particles.pop(r)
    return children, Particles

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
    Simulates the dying of one species, i.e. the reaction  A-> zero. However only
    species in the list 'NotImmune' can die.
    It returns a new list of arrays (position of the species) with removed
    dead species.
    Particles=list of 2D arrays (positions)
    rate=microscopic reaction rate
    deltat=time-step size
    NotImmune= positition of Particles that actually are able to die
    If all particles can die, just set NotImmune=Particles
    '''
    
    pdying=1-np.exp(-rate*deltat) 
    for r in reversed(range(len(Particles))):
        if pdying > np.random.rand() and any((Particles[r] == x).all() for x in NotImmune):
            Particles.pop(r)
    return Particles


def second_order_reaction(A, B, rate, deltat, sigma):
    
    '''Simulates the second order reaction A+B->2B, if A and B are closer then sigma 
    (reaction radius). It returns the new list of A particles with removed particles, a list of 
    B particles with new particles and a list of the new particles (new list of B=previous list of B +children)
    A, B=list of 2D arrays
    rate=microscopic reaction rate
    deltat=time-step size
    sigma=reaction radi
    '''
    p=1-exp(-deltat*rate)
    children=[]
    for i in reversed(range(len(A))):
        for j in reversed(range(len(B))):
            if euclid(A[i], B[j])<= sigma and p>np.random.rand():
                children.append(A[i])
                B.pop(j)
                A.pop(i)
                break
                 
    return A, children, B

def virtual(L, deltar, N,i):
    
    '''
    Assigns to virtual particles at the boundary in cells a position,
    such that they can react with each other in the function 'eatcompact'. Returns a list of 2D arrays.
    L=x-boundary coordinate
    deltar=length of the boundary
    N=number of particles we assign a position to
    i=boundary cell
    '''
    
    Virtual=[None]*N
    for j in range(N):
        Virtual[j]=np.array([np.random.uniform(L-deltar, L), np.random.uniform(deltar*i,deltar*(i+1))])
    return Virtual



def eatcompact(A, B, L, deltar1, deltar2, BC1, BC2, rate, sigma, deltat, Lx):
    '''
    Hybrid algorithm for second order reaction S+I->2I, if S and I are closer than sigma 
    (reaction radius). It returns the new list of S particles with removed particles, a list of 
    I particles with new particles and a list of the new particles (new list of B = previous list of B + children).
    A, B = list of 2D arrays
    deltar = length of boundary domain
    BC1 = boundary concentration of A particles (Prey)
    BC2 = boundary concentration of B particles (Pred)
    rate = microscopic reaction rate
    sigma = reaction radius
    L = x-boundary coordinate
    deltat = time-step size
    Lx = x-boundary coordinate
    '''
    
    kindergarten = []

    # Create all Virtual Predators and Preys
    VirtualPreys = []
    for i in range(len(BC1)): 
        integer = int(BC1[i])
        
        VirtualPreys.extend(virtual(L, deltar1, integer, i))

    VirtualPreds = []
    for i in range(len(BC2)):
        integer = int(BC2[i])
        VirtualPreds.extend(virtual(L, deltar2, integer, i))

    # Combine actual and virtual particles for reactions
    AllPreys = A + VirtualPreys
    AllPreds = B + VirtualPreds

    # Perform reactions
    SurvivedPreys, children, NotReactedPred = second_order_reaction(AllPreys, AllPreds, rate, deltat, sigma)
    
    DecPreys = np.zeros(len(BC1))
    SingleVirtualPrey = []
    for i in range(len(BC1)):
         dec = BC1[i] - int(BC1[i])
         SingleVirtualPrey.append(virtual(L, deltar1, 1, i))
         DecPreys[i] = dec

    DecPreds = np.zeros(len(BC2))
    SingleVirtualPred = []
    for i in range(len(BC2)):
         dec = BC2[i] - int(BC2[i])
         SingleVirtualPred.append(virtual(L, deltar2, 1, i))
         DecPreds[i] = dec

    for i in range(len(SingleVirtualPrey)):
        Ignore, children_virtual, NotReactedPred = second_order_reaction(SingleVirtualPrey[i], NotReactedPred, rate * DecPreys[i], deltat, sigma)
        children.extend(children_virtual)

    for i in range(len(SingleVirtualPred)):
        SurvivedPreys, children_virtual, Ignore = second_order_reaction(SurvivedPreys, SingleVirtualPred[i], rate * DecPreds[i], deltat, sigma)
        children.extend(children_virtual)
    
    # Filter out the virtual particles
    FinalPreys = [p for p in SurvivedPreys if not any(np.array_equal(p, v) for v in VirtualPreys)]
    FinalPreds = [p for p in NotReactedPred if not any(np.array_equal(p, v) for v in VirtualPreds)]
 
    
    return FinalPreys, children, FinalPreds