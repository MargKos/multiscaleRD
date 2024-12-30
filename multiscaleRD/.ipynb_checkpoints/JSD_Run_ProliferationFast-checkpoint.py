from __future__ import division
import numpy as np
from Parameters_Proliferation import *
import numpy as np
from scipy.stats import entropy


# calculates the JSD for B example, Sus refers to Prey
ReferenceSus = np.load('./Solutions/FDSolution_Proliferation.npy')

ts=0.1
h=a/(l-1)

constant=int(ts/deltat)

sim=1200
sim_values = [50,100,200,300,400, 500, 600,700,800,900,1000, 1100,1200]
Times = [39, 69, 89]
batches = 50
xbins = np.arange(0, a + h, h)
ybins = np.arange(0, a + h, h)


def HybridPlot(Average, Reference, bd):
    """
    Creates hybrid plots, where the right side corresponds to the FD solution
    """
    Average_t = np.transpose(Average)  # Get reservoir
    
    # Create matrix from particle discretization
    Hybrid = np.zeros(shape=(l_coupling, l_coupling))
    
    for i in range(l_coupling):
        for j in range(int(l_coupling / 2)):
            Hybrid[i, j] = Average_t[i, j]
        
        for j in range(int(l_coupling / 2)):
            Hybrid[i, j + int(l_coupling / 2)] = Reference[i, j + int(l_coupling / 2)]
    return Hybrid


def Discretization(Particles):
    """
    Return the concentration of particles.
    Particles: list of 2D arrays
    """
    xPositions = [p[0] for p in Particles]
    yPositions = [p[1] for p in Particles]
    concentration, _, _ = np.histogram2d(
        xPositions, yPositions, bins=(xbins, ybins),
        weights=np.ones_like(xPositions) / (h**2)
    )
    return concentration


def functionAverage(Indices, time, all_sus):
    """
    Returns the mean-field concentration for each time-step by averaging over all
    simulations
    """
    num_bins = int(a / h)
    DiscreteSus = np.empty([len(Indices), num_bins, num_bins])
   
    for k, idx in enumerate(Indices):
        # Get preloaded particle data
        discrete_sus = Discretization(all_sus[idx])
       
        DiscreteSus[k, :, :] = discrete_sus
      
    # Return averages over all selected simulations
    return DiscreteSus.mean(axis=0)


def JSD(P, Q):
    """
    Calculate the Jensen-Shannon Divergence between P and Q
    """
    P_t = P.flatten()
    Q_t = Q.flatten()

    # Normalize P and Q
    norm_Q = np.sum(Q)
    P_t = P_t / norm_Q
    Q_t = Q_t / norm_Q

    M = 0.5 * (P_t + Q_t)
    return 0.5 * (entropy(P_t, M) + entropy(Q_t, M))


# Preload all particle data into memory
all_sus_files = [
    np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/ExplicitProliferationParticles_{s}_time{time}.npy', allow_pickle=True)
    for s in range(sim)
    for time in Times
]
     
# Results arrays
SusJSDSimulations = np.zeros((len(sim_values), len(Times)))

# Main loop over fixed `s` values
for s_idx, s in enumerate(sim_values):
    for t_idx, time in enumerate(Times):
        batch_jsd_sus = 0
       

        for _ in range(batches):
            # Randomly sample `s` simulations
            SimulationIndices = np.random.choice(range(sim), size=s, replace=False)

            # Compute average for the sampled simulations
            DiscreteSusAverage = functionAverage(SimulationIndices, time, all_sus_files)

            # Hybrid plots
            HybridSus = HybridPlot(DiscreteSusAverage, ReferenceSus[int(constant * (time + 1))], l_coupling)
            
            # Calculate JSD for each
            batch_jsd_sus += JSD(HybridSus, ReferenceSus[int(constant * (time + 1))])
            
        # Average over batches
        SusJSDSimulations[s_idx, t_idx] = batch_jsd_sus / batches
       
# Save results

np.save('./Solutions/JSD_Proliferation.npy', SusJSDSimulations)

  