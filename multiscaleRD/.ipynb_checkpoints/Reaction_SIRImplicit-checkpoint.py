# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:42:00 2020

@author: Margarita
"""

from __future__ import division
import numpy as np
from math import exp
import time

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
    A=Sus
    B=Inf
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


def eatcompact(A, B, Lx, deltar1, deltar2,BC1, BC2,rate, sigma, deltat,Mdiscr):
    
    # Perform reactions inside the particle domain
    
    AllPreys = A + []
    AllPreds = B + []

    # Perform reactions
    
    kindergarten=[]
    
    
    # Perform reactions across borders A (Sus) in PD and B (Inf) in concentration with the implicit algorithm
    for i in range(len(BC2)):
        
        for j in reversed(range(len(AllPreys))):
            
            # calculate the distance from A to the reservoir
            
            d=abs(AllPreys[j][0]-Lx)
            if d<=sigma:
                h=sigma-d
                Vcap = sigma**2 * np.arccos((sigma - h) / sigma) - (sigma - h) * np.sqrt(2 * sigma * h - h**2)
                
                Lambda1=Vcap*rate*BC2[i]/((deltar2)**2)*(deltat/Mdiscr)
                p=1-np.exp(-Lambda1)
               
            else:
                p=0
            for k in range(Mdiscr):
                if p>np.random.uniform():
                    # A reacts with B in the concentration domain
                    # remove A 
                    
                    AllPreys.pop(j)
                    
                   
                    break  # Break out of the `k` loop
        
         
                
    
    # Perform reactions across borders A (Sus) in concentration and B  (Inf)in particle
    children2 = []
    
    for i in range(len(BC1)):
        
        
        
        for j in reversed(range(len(AllPreds))):
            
            # Calculate the distance from A to the reservoir
            d = abs(AllPreds[j][0] - Lx)
            
            if d <= sigma:
                h = sigma - d
                Vcap = sigma**2 * np.arccos((sigma - h) / sigma) - (sigma - h) * np.sqrt(2 * sigma * h - h**2)
                
                Lambda = Vcap * rate * BC1[i] / ((deltar1)**2) * (deltat / Mdiscr) 
                
                p = 1 - np.exp(-Lambda)
            else:
                p = 0
            
            for k in range(Mdiscr):
                
                if p > np.random.uniform():
                    # A reacts with B in the concentration domain
                    # Remove A
                    children2.append(AllPreds[j])
                    AllPreds.pop(j)
                    
                    
                    break  # Break out of the `k` loop
    
          

    
    # reactions in PB domain
    kindergarten.extend(children2)            
    # Pause the timer for internal reactions
    internal_start = time.time()
    AllPreys, childrenIntern, AllPredsAfterIntern = second_order_reaction(AllPreys, AllPreds, rate, deltat, sigma, Lx)
    kindergarten.extend(childrenIntern)
    internal_end = time.time()

    # Calculate elapsed time
    internal_time = internal_end - internal_start
 
   

    return AllPreys, kindergarten,AllPredsAfterIntern, internal_time



