from __future__ import division
import numpy as np
from Parameters_LV import *
import numpy as np
from scipy.stats import entropy


ReferencePrey = np.load('./Solutions/FDLVSolution1.npy')  # gets Data from continuous solution of A
ReferencePred  = np.load('./Solutions/FDLVSolution2.npy')  # gets Data from continuous solution of B

# Parameters
ts = 0.1
constant = int(ts / deltat)
sim_values = [1,10,25,50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
Times = [9, 14,24]
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


def functionAverage(Indices, time, all_sus, all_inf):
    """
    Returns the mean-field concentration for each time-step by averaging over all
    simulations
    """
    num_bins = int(a / h)
    DiscreteSus = np.empty([len(Indices), num_bins, num_bins])
    DiscreteInf = np.empty([len(Indices), num_bins, num_bins])
   

    for k, idx in enumerate(Indices):
        # Get preloaded particle data
        discrete_sus = Discretization(all_sus[idx])
        discrete_inf = Discretization(all_inf[idx])
        
        DiscreteSus[k, :, :] = discrete_sus
        DiscreteInf[k, :, :] = discrete_inf
        
    # Return averages over all selected simulations
    return DiscreteSus.mean(axis=0), DiscreteInf.mean(axis=0)


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
all_prey_files = [
    np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/LVPreyParticles{s}time{time}.npy', allow_pickle=True)
    for s in range(500)
    for time in Times
]
all_pred_files = [
    np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/LVPredatorParticles{s}time{time}.npy', allow_pickle=True)
    for s in range(500)
    for time in Times
]        
# Results arrays
PreyJSDSimulations = np.zeros((len(sim_values), len(Times)))
PredJSDSimulations = np.zeros((len(sim_values), len(Times)))

# Main loop over fixed `s` values
for s_idx, s in enumerate(sim_values):
    for t_idx, time in enumerate(Times):
        batch_jsd_Prey = 0
        batch_jsd_Pred = 0

        for _ in range(batches):
            # Randomly sample `s` simulations
            SimulationIndices = np.random.choice(range(500), size=s, replace=False)
            print(SimulationIndices, s)

            # Compute average for the sampled simulations
            DiscretePreyAverage, DiscretePredAverage = functionAverage(
                SimulationIndices, time, all_prey_files, all_pred_files
            )

            # Hybrid plots
            HybridPrey = HybridPlot(DiscretePreyAverage, ReferencePrey[int(constant * (time + 1))], l_coupling)
            HybridPred = HybridPlot(DiscretePredAverage, ReferencePred[int(constant * (time + 1))], l_coupling)
           
            # Calculate JSD for each
            batch_jsd_Prey += JSD(HybridPrey, ReferencePrey[int(constant * (time + 1))])
            batch_jsd_Pred += JSD(HybridPred, ReferencePred[int(constant * (time + 1))])
           

        # Average over batches
        PreyJSDSimulations[s_idx, t_idx] = batch_jsd_Prey / batches
        PredJSDSimulations[s_idx, t_idx] = batch_jsd_Pred / batches

#%%
np.save('./Solutions/LVPreyJSDPaper.npy', PreyJSDSimulations)
np.save('./Solutions/LVPredJSDPaper.npy', PredJSDSimulations)
        

