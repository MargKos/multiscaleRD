from __future__ import division
import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


'B Mauricio Parameters Number 2 (not in the paper)'

def function1(x):
	if 5.1<=x[0]<=6.1 and 1.5<=x[1]<=2.5:
		ux=1200 
	elif 5.1<=x[0]<=6.1 and 4.5<=x[1]<=5.5:
		ux=1200 
	elif 5.1<=x[0]<=6.1 and 7.5<=x[1]<=8.5:
		ux=1200 
	else:
		ux=0
	
	return ux    

def function2(x):
    
    if np.linalg.norm(x-np.array([6.5,2.5]))<0.5: 
        ux=800 

    elif np.linalg.norm(x-np.array([6.5,7.5]))<0.5: 
        ux=800 
       
    else:
        ux=0
    
    return ux


l=100+1
m=l
deltat=0.005
timesteps=300+1
a=10
h=a/(l-1) # grid size 

maxtime=(timesteps-1)*deltat

D1=0.8 
D2=0.8 
D3=0.8

sigma=0.01
r1_macro=0.015 # infection rate macro
r2=0.2 # recovery rate
r2_macro=r2
r1=r1_macro/(np.pi*(sigma**2))  # infection rate mictoe
Lx=a/2 # half of the grid
deltar1=np.sqrt(deltat*D1*2) # boundary size
deltar2=np.sqrt(deltat*D2*2) # boundary size
deltar3=np.sqrt(deltat*D3*2) # boundary size

l_coupling=l-1 
save_time=0.1 