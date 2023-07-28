import numpy as np
from scipy.spatial import Voronoi
import selfpropelledparticlevoronoi as sppv
import matplotlib.pyplot as plt
from scipy import ndimage as nd 
import findconnections as fc

def newMod(a,b):
    res = np.zeros(a.shape)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            res[i,j] = a[i,j]%b
            if a[i,j]<0:
                res[i,j] -= b
    return res 

def newWhere(a,b):
    newA = np.zeros(a.shape)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            newA[i,j] = np.abs(a[i,j])<b
            if np.abs(a[i,j])>b:
                newA[i,j] =  - 1
    return newA 


def thermalize(regions, point_region,cells,vel,r0):
    LC = len(cells)
    FC = []
    
    for c in range(LC):
        xy_c = cells[c]
        neigh_c = fc.find_center_neighbour_center(regions,point_region,c)
        f_c = 0.0*np.array([1.,1.])#+(np.random.rand(2)-0.5)*2
        for n_c in neigh_c:
            xy_n = cells[n_c]
            v_nc = xy_c-xy_n
            r_nc = fc.norm(v_nc)
            l_nc =v_nc/(r_nc+1e-1)
            
            #Lennard Jones
            
            rho = r0/(r_nc+1e-1)
            f_c += -3*(rho**2-rho)*l_nc-(r_nc-1.25*r0)*l_nc
                
        FC.append(f_c)         
        
    return np.array(FC)


def newcoords(N):
    L_max = 4.5
    L_min = -4.5

    x_coord = (L_max-L_min)*(np.random.rand(N)) + L_min
    y_coord = (L_max-L_min)*(np.random.rand(N)) + L_min

    coords = np.array((x_coord,y_coord)).T

    vel_array = np.zeros((coords.shape))

    #Prior to simulation run a thermalization process

    thresh_f = 1e-1

    avg_f = 1

    dt = 5*1e-2
    DeltaL = L_max - L_min
    r0 = min(DeltaL/(3*np.sqrt(N)),0.5)

    steps = 0
    
    n_steps = int(input('Max steps for thermalization: '))
    
    fft_array = np.zeros((100,n_steps))


    while (avg_f > thresh_f) and steps < n_steps:
        
        vor = Voronoi(coords)
        vorPointRegion = vor.point_region
        vorRegions = vor.regions
    
        F_center = thermalize(vorRegions,vorPointRegion,coords,vel_array,r0)
        vel_array = F_center
    
        #Periodic boundary conditions
    
        #coords = newMod(coords + vel_array*dt,5)
    
        #Reflexive boundary conditions
    
        A = newWhere(coords + vel_array*dt,5)
        coords = coords + vel_array*dt*A

        avg_f = np.mean(F_center**2)**0.5
        
        steps += 1
        
    print(steps)
    print(avg_f)
    print("done")
    return coords

