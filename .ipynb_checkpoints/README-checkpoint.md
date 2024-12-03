# MultiscaleRD

This repository contains code to accompany the paper *"Coupling particle-based reaction-diffusion simulations with reservoirs mediated by reaction-diffusion PDEs"* by M. Kostré, C. Schütte, F. Noé, and M. J. del Razo.  
[[arXiv]](https://arxiv.org/pdf/2006.00003.pdf)

In this code, we implement a hybrid scheme that couples particle-based reaction-diffusion simulations with spatially and time-dependent reservoirs mediated by reaction-diffusion PDEs. It solves the reaction-diffusion PDE using a finite difference scheme with a Crank-Nicolson integrator, and implements particle-based reaction-diffusion simulations based on the Doi model, similar to how it is done in [ReaDDy2](https://readdy.github.io/).  
The hybrid scheme consistently couples the particle-based simulation to the PDE-mediated reservoir. We verify the scheme using **three examples**:

1. **Diffusion process in 1D**: .
2. **Diffusion process with a proliferation reaction**: This can be generalized to systems with zeroth- and/or first-order reactions.  (`multiscaleRD/FD_Proliferation.py` and `multiscaleRD/Coupling_Proliferation.py`)
3. **Lotka-Volterra (predator-prey) reaction-diffusion process**: This can be generalized to systems with up to second-order reactions. (`multiscaleRD/FD_LV.py` and `multiscaleRD/Coupling_LV.py`)
4. **SIR (susceptible, infected, recovered) reaction-diffusion process**: This can be generalized to systems with up to second-order reactions. (`multiscaleRD/FD_SIR.py` and `multiscaleRD/Coupling_SIR.py`)

## Requirements

- Python 3.x
- Numpy 1.24.2
- Matplotlib 3.63 or above
- Pickle (part of Python standard library)
- Git (to clone this repository)

## How to run this code?

1. Clone this repository to your local machine:  
   `git clone https://github.com/MargKos/multiscaleRD.git`
2. (Optional) Modify the default mathematical and numerical parameters for the finite difference scheme and particle-based simulation in `multiscaleRD/Parameters_(example).py`.
3. Solve the reaction-diffusion PDE for each example by running the finite difference code (`multiscaleRD/FD_(example).py`). This will generate the reference solution(s) and reservoir.
4. Run multiple particle-based simulations by executing `multiscaleRD/Coupling_(example).py`, preferably across multiple computers. For each time-step, store the locations of particles from each simulation.
5. Use `multiscaleRD/Discretization_(example).py` to calculate the average over the particle-based simulations.
6. Run `multiscaleRD/Plot_(example).py` to plot the concentration.
7. Use `multiscaleRD/CompareMean_(example).py` to compare the number of particles across the particle-based domain.
8. To calculate the Jensen-Shannon Divergence (JSD), run `multiscaleRD/JSD_Run_(example).py` and plot it using `multiscaleRD/JSD_Plot_(example).py`.

## Folder Organization

- `multiscaleRD/Parameters_(example).py`: Contains mathematical parameters (reaction rates, diffusion coefficients), numerical parameters (time-step size, grid size, boundary size), and computational parameters (number of parallel simulations).
- `multiscaleRD/FD_(example).py`: Implements the finite difference scheme for the reaction-diffusion equation with homogeneous Neumann boundary conditions.
- `multiscaleRD/Coupling_(example).py`: Contains the hybrid algorithm.
- `multiscaleRD/Injection.py` and `multiscaleRD/Reaction.py`: Contain functions for the injection and reaction procedures used in `multiscaleRD/Coupling_(example).py`:
  - `multiscaleRD/Reaction_(example).py` implements particle-based reaction-diffusion simulations (up to second-order reactions), such as proliferation (A → 2A), degradation (A → 0), and diffusion of particles modeled by the Euler-Maruyama scheme for a) LV and proliferation examples b) SIR example. These functions return, for example, lists of newly created particle positions.
- The folder `multiscaleRD/Simulations`: Contains all the particle trajectories resulting from the hybrid scheme simulations (`multiscaleRD/Coupling_(example).py`).
The second and third example are using Reaction_LV and Injection_LV functions, where SIR uses Reaction_SIR and Injection_SIR functions
- The folder `multiscaleRD/Solutions`: Contains solutions of the reaction-diffusion PDE and the average over particle-based simulations from `multiscaleRD/Simulations`. These files are used for plotting or analysis in `multiscaleRD/Plot_(example).py`.

## Notes on Code Notation

- The reference solutions of the preys/susceptible species are denoted by "1", the solutions of the predators/infected by "2", and the recovered by "3".

## Sample Solutions

- LV
This example shows two videos for the reaction-diffusion dynamics of the concentration of prey in the Lotka-Volterra system (predator-prey). The left video shows the reference simulation obtained with a finite-difference scheme. The right video shows the results of the hybrid simulation. Here, the left half corresponds to the average concentration over several particle-based simulations using our scheme, and the right half corresponds to the PDE-mediated reservoir.

<img src="Videos/ReferencePrey_video.gif" width="400"> <img src="Videos/HybridPrey_video.gif" width="400" />

<img src="Videos/ReferencePred_video.gif" width="400"> <img src="Videos/HybridPred_video.gif" width="400" />

- SIR


This example shows two videos for the SIR dynamics of the concentration of 1) susceptible 2) infected and 3) recovered species in the system for setting B. The left video shows the reference simulation obtained with a finite-difference scheme. The right video shows the results of the hybrid simulation with the tau-leaping scheme. Here, the left half corresponds to the average concentration over several particle-based simulations using our scheme, and the right half corresponds to the PDE-mediated reservoir.

<img src="Videos/ReferenceSus_video_B.gif" width="400"> <img src="Videos/HybridSus_video_TauB.gif" width="400" />

<img src="Videos/ReferencInf_video_B.gif" width="400"> <img src="Videos/HybridInf_video_TauB.gif" width="400" />

<img src="Videos/ReferenceRec_video_B.gif" width="400"> <img src="Videos/HybridRec_video_TauB.gif" width="400" />

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.