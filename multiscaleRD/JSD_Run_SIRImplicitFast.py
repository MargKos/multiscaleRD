from __future__ import division
import numpy as np
from ParametersB_SIR import *
from scipy.stats import entropy

# Preload reference data into memory
ReferenceSus = np.load('./Solutions/FDSIR1_B.npy')
ReferenceInf = np.load('./Solutions/FDSIR2_B.npy')
ReferenceRec = np.load('./Solutions/FDSIR3_B.npy')

# Parameters
ts = 0.1
constant = int(ts / deltat)
sim_values = [1,10,25,50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
Times = [4, 9]
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


def functionAverage(Indices, time, all_sus, all_inf, all_rec):
    """
    Returns the mean-field concentration for each time-step by averaging over all
    simulations
    """
    num_bins = int(a / h)
    DiscreteSus = np.empty([len(Indices), num_bins, num_bins])
    DiscreteInf = np.empty([len(Indices), num_bins, num_bins])
    DiscreteRec = np.empty([len(Indices), num_bins, num_bins])

    for k, idx in enumerate(Indices):
        # Get preloaded particle data
        discrete_sus = Discretization(all_sus[idx])
        discrete_inf = Discretization(all_inf[idx])
        discrete_rec = Discretization(all_rec[idx])

        DiscreteSus[k, :, :] = discrete_sus
        DiscreteInf[k, :, :] = discrete_inf
        DiscreteRec[k, :, :] = discrete_rec

    # Return averages over all selected simulations
    return DiscreteSus.mean(axis=0), DiscreteInf.mean(axis=0), DiscreteRec.mean(axis=0)


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
    np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIRSusImplicitB{s}time{time}.npy', allow_pickle=True)
    for s in range(500)
    for time in Times
]
all_inf_files = [
    np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIRInfImplicitB{s}time{time}.npy', allow_pickle=True)
    for s in range(500)
    for time in Times
]
all_rec_files = [
    np.load(f'/home/htc/bzfkostr/SCRATCH/SimulationsMultiscale/SIRRecoveryImplicitB{s}time{time}.npy', allow_pickle=True)
    for s in range(500)
    for time in Times
]

# Results arrays
SusJSDSimulations = np.zeros((len(sim_values), len(Times)))
InfJSDSimulations = np.zeros((len(sim_values), len(Times)))
RecJSDSimulations = np.zeros((len(sim_values), len(Times)))

# Main loop over fixed `s` values
for s_idx, s in enumerate(sim_values):
    for t_idx, time in enumerate(Times):
        batch_jsd_sus = 0
        batch_jsd_inf = 0
        batch_jsd_rec = 0

        for _ in range(batches):
            # Randomly sample `s` simulations
            SimulationIndices = np.random.choice(range(500), size=s, replace=False)

            # Compute average for the sampled simulations
            DiscreteSusAverage, DiscreteInfAverage, DiscreteRecAverage = functionAverage(
                SimulationIndices, time, all_sus_files, all_inf_files, all_rec_files
            )

            # Hybrid plots
            HybridSus = HybridPlot(DiscreteSusAverage, ReferenceSus[int(constant * (time + 1))], l_coupling)
            HybridInf = HybridPlot(DiscreteInfAverage, ReferenceInf[int(constant * (time + 1))], l_coupling)
            HybridRec = HybridPlot(DiscreteRecAverage, ReferenceRec[int(constant * (time + 1))], l_coupling)

            # Calculate JSD for each
            batch_jsd_sus += JSD(HybridSus, ReferenceSus[int(constant * (time + 1))])
            batch_jsd_inf += JSD(HybridInf, ReferenceInf[int(constant * (time + 1))])
            batch_jsd_rec += JSD(HybridRec, ReferenceRec[int(constant * (time + 1))])

        # Average over batches
        SusJSDSimulations[s_idx, t_idx] = batch_jsd_sus / batches
        InfJSDSimulations[s_idx, t_idx] = batch_jsd_inf / batches
        RecJSDSimulations[s_idx, t_idx] = batch_jsd_rec / batches

# Save results
np.save('./Solutions/JSDSusImplicitB.npy', SusJSDSimulations)
np.save('./Solutions/JSDRecImplicitB.npy', RecJSDSimulations)
np.save('./Solutions/JSDInfImplicitB.npy', InfJSDSimulations)
