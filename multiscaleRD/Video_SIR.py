#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 12:11:01 2019

@author: bzfkostr
"""

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from ParametersA_SIR import *
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.animation import FuncAnimation

ts=0.1 # discretization length of results
constant=int(ts/deltat)
DiscreteSusAverage = np.load('./Solutions/SusAverageB.npy')
DiscreteInfAverage = np.load('./Solutions/InfAverageB.npy')
DiscreteRecAverage = np.load('./Solutions/RecAverageB.npy')

ReferenceSus = np.load('./Solutions/FDSIR1_B.npy')
ReferenceInf = np.load('./Solutions/FDSIR2_B.npy')
ReferenceRec = np.load('./Solutions/FDSIR3_B.npy')


timesteps_cut=int(deltat*(timesteps-1)/ts)

def HybridPlot(Average, Reference, bd):
    '''
    Creates hybrid plots, where the right side corresponds to the FD solution
    and the left side to the mean-field concentration obtained from the coupling.
    Average=mean-field concentration
    Concentration=FD solution
    bd=location of the boundary
    '''

    listH = []  # list of Hybrid solutions
    listR = []  # list of Reference solutions
    
    constant=int(ts/deltat)
    Hmean=np.zeros(timesteps_cut)
    Cmean=np.zeros(timesteps_cut)
    if len(Average)!=timesteps_cut:
        print('error')
        
    for t in range(timesteps_cut):
        Average_t = np.transpose(Average[t])  # get reservoir
        Ref_t=Reference[(t+1)*constant]
        #Ref_t=Reference[(t)*constant]
        
        # create matrix fromt he particle discretization
        Particle = np.zeros(shape=(l_coupling,l_coupling))
        
        for i in range(l_coupling):
            for j in range(int(l_coupling / 2) ):
                Particle[i, j] = Average_t[i, j]
            
            for j in range(int(l_coupling / 2)):
                
                Particle[i,j+int(l_coupling/2)]=Ref_t[i,j+int(l_coupling/2)]
        listH.append(Particle)
            
        # average
        Hmean[t]=np.mean(Particle)
            
         
    for t in range(timesteps_cut):

        Ref_t=Reference[(t+1)*constant]
        #for t in range(timesteps_cut):
        #Ref_t=Reference[t*constant]

        Cmean[t]=np.mean(Ref_t)
        listR.append(Ref_t)

    return Hmean, Cmean, listH, listR


sus_mean, RefS, HybridSus, ReferenceSus = HybridPlot(DiscreteSusAverage, ReferenceSus, l_coupling)
inf_mean, RefI,  HybridInf, ReferenceInf= HybridPlot(DiscreteInfAverage,  ReferenceInf, l_coupling)
rec_mean, RefR,  HybridRec, ReferenceRec = HybridPlot(DiscreteRecAverage, ReferenceRec, l_coupling)


# Define the custom colormap with specified colors
colors1 = ['black','indigo', 'royalblue', 'lightskyblue', 'azure']
n_bins1 = 256  # Number of bins for the colormap

# Create the colormap
custom_colormap1 = LinearSegmentedColormap.from_list('custom_colormap', colors1, N=n_bins1)



# Updated function for creating animation
def init_plot(Data, Max):
    '''Initializes the plot for animation.'''
    fig, ax = plt.subplots(figsize=(a, a))
    im = ax.imshow(Data[0], interpolation='nearest', cmap=custom_colormap1, extent=[0, a, 0, a])
    ax.set_xlabel('x', fontsize=26)
    ax.set_ylabel('y', fontsize=26)
    ax.tick_params(labelsize=26)
    cbar = plt.colorbar(im, ax=ax, fraction=0.045)
    cbar.ax.tick_params(labelsize=26)
    im.set_clim(-Max / 20, Max)
    plt.tight_layout()
    return fig, ax, im

# Update function for each frame in the animation
def update_frame(frame, Data, im):
    '''Updates the frame of the plot for animation.'''
    im.set_array(Data[frame])
    return [im]


MaxPlot=[200,400,300]

# Generate video for Hybrid Susceptible
fig, ax, im = init_plot(HybridSus, MaxPlot[0])
anim = FuncAnimation(fig, update_frame, frames=len(HybridSus), fargs=(HybridSus, im), blit=True)

# Save animation as a video file
anim.save('./Videos/HybridSus_video_B.gif', writer='ffmpeg', fps=10)

# Generate video for Reference Susceptible
fig, ax, im = init_plot(ReferenceSus, MaxPlot[0])
anim_ref = FuncAnimation(fig, update_frame, frames=len(ReferenceSus), fargs=(ReferenceSus, im), blit=True)

# Save animation as a video file
anim_ref.save('./Videos/ReferenceSus_video_B.gif', writer='ffmpeg', fps=10)


# Generate video for Hybrid Infected
fig, ax, im = init_plot(HybridInf, MaxPlot[1])
anim = FuncAnimation(fig, update_frame, frames=len(HybridInf), fargs=(HybridInf, im), blit=True)

# Save animation as a video file
anim.save('./Videos/HybridInf_video_B.gif', writer='ffmpeg', fps=10)

# Generate video for Reference Infected
fig, ax, im = init_plot(ReferenceInf, MaxPlot[1])
anim_ref = FuncAnimation(fig, update_frame, frames=len(ReferenceInf), fargs=(ReferenceInf, im), blit=True)

# Save animation as a video file
anim_ref.save('./Videos/ReferenceInf_video_B.gif', writer='ffmpeg', fps=10)

# Generate video for Hybrid Infected
fig, ax, im = init_plot(HybridRec, MaxPlot[2])
anim = FuncAnimation(fig, update_frame, frames=len(HybridRec), fargs=(HybridRec, im), blit=True)

# Save animation as a video file
anim.save('./Videos/HybridRec_video_B.gif', writer='ffmpeg', fps=10)

# Generate video for Reference Infected
fig, ax, im = init_plot(ReferenceRec, MaxPlot[2])
anim_ref = FuncAnimation(fig, update_frame, frames=len(ReferenceRec), fargs=(ReferenceRec, im), blit=True)

# Save animation as a video file
anim_ref.save('./Videos/ReferenceRev_video_B.gif', writer='ffmpeg', fps=10)

