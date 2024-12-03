from __future__ import division
import numpy as np
from Parameters_LV2 import *
import numpy as np
from scipy.stats import entropy


ReferencePrey = np.load('./Solutions/PaperFDLVSolution1.npy')  # gets Data from continuous solution of A
ReferencePred  = np.load('./Solutions/PaperFDLVSolution2.npy')  # gets Data from continuous solution of B

ts=0.1
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

    DiscretePred=np.empty([sim, number_bins,number_bins])
    DiscretePredAverage=np.empty([ number_bins,number_bins])
    
    k=0
    for s in Indices:

        Prey_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/'+str('PaperLVPreyParticles')+str(s)+ 'time'
                +str(time)+'.npy', allow_pickle=True)
        Predator_s_t=np.load('/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/'+str('PaperLVPredatorParticles')+str(s)+ 'time'
                +str(time)+'.npy', allow_pickle=True)
       

        discrete_Prey_s_t=Discretization(a, l_coupling, Prey_s_t)
        
        DiscretePrey[k,:,:]=discrete_Prey_s_t

        # same for predator

        discrete_Predator_s_t=Discretization(a, l_coupling, Predator_s_t)

        DiscretePred[k,:,:]=discrete_Predator_s_t
        
        k=k+1
            

    DiscretePreyAverage=np.mean(DiscretePrey[:,:,:], axis=0)
    DiscretePredAverage=np.mean(DiscretePred[:,:,:], axis=0)
  
        
    return DiscretePreyAverage, DiscretePredAverage



sim=500
Times=[39,69,89]

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
PredJSDSimulations=np.zeros((sim, len(Times)))

batches=50
for t in range(len(Times)):
    for s in range(sim):
        for j in range(batches):
            SimulationIndices=np.random.randint(0,sim, s) # pick s simulations in the range from 0,500

            DiscretePreyAverage, DiscretePredAverage=functionAverage(SimulationIndices, Times[t])
                
            HybridPrey = HybridPlot(DiscretePreyAverage, ReferencePrey[int(constant*(Times[t]+1))], l_coupling)
            HybridPred= HybridPlot(DiscretePredAverage,  ReferencePred[int(constant*(Times[t]+1))], l_coupling)
          
            P=HybridPrey
            Q=ReferencePrey[int(constant*(Times[t]+1))]
            jsd=JSD(P, Q)
            PreyJSDSimulations[s,t]=jsd+PreyJSDSimulations[s,t]        
            
            P=HybridPred
            Q=ReferencePred[int(constant*(Times[t]+1))]
            jsd=JSD(P, Q)
            PredJSDSimulations[s,t]=jsd+PredJSDSimulations[s,t]        

    
        PreyJSDSimulations[s,t]=PreyJSDSimulations[s,t]/(batches+1)
        PredJSDSimulations[s,t]=PredJSDSimulations[s,t]/(batches+1)
                                                 
                                                 
    
#%%
np.save('./Solutions/LVPreyJSDPaper.npy', PreyJSDSimulations)
np.save('./Solutions/LVPredJSDPaper.npy', PredJSDSimulations)
        

