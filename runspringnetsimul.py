import numpy as np
from scipy.spatial import Voronoi
import selfpropelledparticlevoronoi as sppv
import pathlib
import h5py as h5
import matplotlib.pyplot as plt
import randomlatticethermalized as rlt

''' 

Run simulations in this module

'''

option_1 = int(input(" Do you want to run the simulation on a random lattice: y(1) or n(0): "))
if option_1 == 0:
    # Define a regular lattice
    N_m = int(input("Number of points you want to add onto the network"))
    
    L_max = 5
    L_min = -5  
    dL = (L_max-L_min)/N_m
    A_0 = np.array([[1,0],[-1/2,np.sqrt(3)/2]])
    A = 1/(2*np.linalg.det(A_0))*(A_0.T+A_0)
    A_inv = np.linalg.inv(A)

    LxyM = A_inv.dot(np.array([L_max,L_max]))
    Lxym = A_inv.dot(np.array([L_min,L_min]))

    u = np.linspace(Lxym[0],LxyM[0],N_m)
    v = np.linspace(Lxym[1],LxyM[1],N_m)

    coord_list = []

    for j in range(N_m):
        for i in range(N_m):
            coord_list.append(A.dot(np.array([u[j],v[i]])))
        
    coords = np.array(coord_list)
    N = N_m**2
    


if option_1 == 1:
    # Define a random lattice

    N = int(input("Square root of the number of points you want to add onto the network "))**2

    coords = rlt.newcoords(N)


# Define arrays to insert numerical results

L_max = 5
L_min = -5  
dt = 5*1e-3

M = 5000

vor = Voronoi(coords)
vorPointRegion = vor.point_region
vorRegions = vor.regions
DeltaL = L_max - L_min
r0 = min(DeltaL/(np.sqrt(N)),0.5)
coords_evo = np.zeros((N,2,M))
coords_evo_vertex = np.zeros((len(vor.vertices),2,M))
coords_evo[:,:,0] = coords
coords_evo_vertex[:,:,0] = vor.vertices

coords_evo_v = np.zeros((N,2,M))
coords_evo_vertex_v = np.zeros((len(vor.vertices),2,M))

vorRidges = sppv.remove_minus(vor.ridge_vertices)

Force_center_vector = np.zeros((N,2,M))
Force_vertex_vector = np.zeros((len(vor.vertices),2,M))

vorRidges_array = np.zeros((len(vorRidges),2,M))
vorRidges_array[:,:,0] = vorRidges

vorRegions_array = []
vorRegions_array.append(vorRegions)

boundary = sppv.find_boundary_vertices(len(vor.vertices),vorRidges)

continue_option = int(input('Continue with simulation: (y-1/n-0): '))

if continue_option == 1:
    
    # Run simulations

    for i in range(M-1):
        #Do T1 transitions and compute forces afterwards
    
        
        #vorRidges_array[:,:,i], coords_evo_vertex[:,:,i]
        F_center, F_vertex, dist_med_v = sppv.force_vtx_elastic(vorRegions, vorPointRegion, vorRidges_array[:,:,i], boundary, 5,r0,coords_evo[:,:,i],coords_evo_vertex[:,:,i])
        
        Force_center_vector[:,:,i] = F_center
        Force_vertex_vector[:,:,i] = F_vertex
        
        coords_evo[:,:,i+1] = coords_evo[:,:,i] + F_center*dt
        coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + F_vertex*dt
        #vorRidges_array[:,:,i+1] = vorRidges_array[:,:,i]
        vorRidges_array[:,:,i+1], coords_evo_vertex[:,:,i+1], vorRegions, transition_counter =  sppv.T1transition(coords_evo_vertex[:,:,i+1],vorRidges_array[:,:,i],vorRegions, vorPointRegion,max(r0/4,dist_med_v/4))
        
        if transition_counter > 0:
            with open('transition_times.txt','a') as text:
                text.write('At time ')
                text.write(str(i))
                text.write(' number of transitions: ')
                text.write(str(transition_counter)+ '\n')
        
        vorRegions_array.append(vorRegions)
            
    main_path = pathlib.Path().absolute()
    datafilesavename = input('Name of savefile: ')
    datafilesavename = datafilesavename + '.h5'

    data_save = h5.File(datafilesavename,'w')
    data_save.create_dataset('centers',data = coords_evo)
    data_save.create_dataset('vertices',data=coords_evo_vertex)
    data_save.create_dataset('forces_center',data=Force_center_vector)
    data_save.create_dataset('forces_vertices',data=Force_vertex_vector)
    data_save.create_dataset('Edge connections',data=vorRidges_array)
    data_save.create_dataset('number_of_points',data=N)
    data_save.create_dataset('time_of_simul',data=M)
    data_save.create_dataset('timestep',data = dt)
    data_save.close()

