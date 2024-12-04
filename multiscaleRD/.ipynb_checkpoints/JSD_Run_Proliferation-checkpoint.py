from __future__ import division
import numpy as np
from Parameters_Proliferation import *
import numpy as np
from scipy.stats import entropy

ReferencePrey = np.load('./Solutions/FDSolution_Proliferation.npy')

timesteps=1000+1
ts=0.1
h=a/(l-1)
timesteps_cut = int(deltat * (timesteps-1) / ts)
print(timesteps_cut, deltat, timesteps-1, ts)
constant=int(ts/deltat)


def HybridPlot(Average, Reference, bd):
    '''
    Creates hybrid plots, where the right side corresponds to the FD solution
    '''
    Average_t = np.transpose(Average)  # get reservoir
        
    # create matrix fromt he particle discretization
    Hybrid = np.zeros(shape=(l_coupling,l_coupling))
        
    for i in range(l_coupling):
        for j in range(int(l_coupling / 2) ):
            Hybrid[i, j] = Average_t[i, j]
        
        for j in range(int(l_coupling / 2)):
            
            Hybrid[i,j+int(l_coupling/2)]=Reference[i,j+int(l_coupling/2)]
    return Hybrid

def Discretization(a, discret, Particles): 
    
    '''
    Return the concentration of particles.
    a=domain length
    discret=discrtization parameter (number of cells)
    Particles=list of 2D arrays
    '''
    
    xPositions=[]
    yPositions=[]
    for i in range(len(Particles)):
       
        xPositions.append(Particles[i][0])
        yPositions.append(Particles[i][1])
    
    xbins=np.arange(0,a+h, h)
    ybins=np.arange(0,a+h, h)
    concentration, xbins, ybins=np.histogram2d(xPositions, yPositions, bins=(xbins, ybins), weights=np.ones_like(xPositions)/(h**2))
    
    return concentration
    

def functionAverage(Indices, time): 
    
    '''
    Returns the mean-field concentration for each time-step by averaging over all
    simulations
    PreySimulation=list of all simulations, see Coupling.py
    '''
    sim=len(Indices)
    
    k=0
    
    
    number_bins=int(a/h)
    
    DiscretePrey=np.empty([sim, number_bins,number_bins])
    DiscretePreyAverage=np.empty([number_bins,number_bins])
    
    k=0
    for s in Indices:    
            
        Prey_s_t=np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/ProliferationParticles_{s}_time{time}.npy', allow_pickle=True)

        discrete_Prey_s_t=Discretization(a, l_coupling, Prey_s_t)
        
        DiscretePrey[k,:,:]=discrete_Prey_s_t

        
        k=k+1

    DiscretePreyAverage=np.mean(DiscretePrey[:,:,:], axis=0)
    
        
    return DiscretePreyAverage



sim=500
Times=[4,9,14]

def JSD(P, Q):

    # go over all time steps
    # Q should be the reference solution
    P_t = P.flatten()
    Q_t = Q.flatten()
    
    # normalize P and Q
    
    norm_Q=np.sum(Q)
    
    P_t=P_t/norm_Q
    Q_t=Q_t/norm_Q

    # calculate JSD

    M = 0.5 * (P_t + Q_t)

    JSD = 0.5 * (entropy(P_t, M) + entropy(Q_t, M))

       

    return JSD

k=0
PreyJSDSimulations=np.zeros((sim, len(Times)))

batches=50
for t in range(len(Times)):
    for s in range(sim):
        for j in range(batches):
            SimulationIndices=np.random.randint(0,sim, s) # pick s simulations in the range from 0,500

            DiscretePreyAverage=functionAverage(SimulationIndices, Times[t])
                
            HybridPrey = HybridPlot(DiscretePreyAverage, ReferencePrey[int(constant*(Times[t]+1))], l_coupling)
           
            P=HybridPrey
            Q=ReferencePrey[int(constant*(Times[t]+1))]
            jsd=JSD(P, Q)
            PreyJSDSimulations[s,t]=jsd+PreyJSDSimulations[s,t]        
   
        PreyJSDSimulations[s,t]=PreyJSDSimulations[s,t]/(batches+1)
                                                 
                                                 
    
#%%
np.save('./Solutions/JSD_Proliferation.npy', PreyJSDSimulations)

        
#%%
