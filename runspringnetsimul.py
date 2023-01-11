import numpy as np
from scipy.spatial import Voronoi
import selfpropelledparticlevoronoi as sppv
import pathlib
import h5py as h5

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

    N = int(input("Square root of the number of points you want to add onto the network"))**2

    L_max = 5
    L_min = -5

    x_coord = (L_max-L_min)*(np.random.rand(N)) + L_min
    y_coord = (L_max-L_min)*(np.random.rand(N)) + L_min

    coords = np.array((x_coord,y_coord)).T

vor = Voronoi(coords)

# Define arrays to insert numerical results

M = 500
coords_evo = np.zeros((N,2,M))
coords_evo_vertex = np.zeros((len(vor.vertices),2,M))
coords_evo[:,:,0] = coords
coords_evo_vertex[:,:,0] = vor.vertices

coords_evo_v = np.zeros((N,2,M))
coords_evo_vertex_v = np.zeros((len(vor.vertices),2,M))

dt = 1.5*1e-1
DeltaL = L_max - L_min
r0 = min(DeltaL/(5*np.sqrt(N)),0.5)

vor_ridges = sppv.remove_minus(vor.ridge_vertices)

Force_center_vector = np.zeros((N,2,M))
Force_vertex_vector = np.zeros((len(vor.vertices),2,M))

vor_ridges_array = np.zeros((len(vor_ridges),2,M))
vor_ridges_array[:,:,0] = vor_ridges

vor_point_region = vor.point_region
vor_regions = vor.regions

vor_regions_array = []
vor_regions_array.append(vor_regions)



# Run simulations

for i in range(M-1):
    
    #Do T1 transitions and compute forces afterwards
 
    
    #vor_ridges_array[:,:,i], coords_evo_vertex[:,:,i]
    F_center, F_vertex, dist_med_v = sppv.force_vor_elastic(vor_regions, vor_point_region, vor_ridges_array[:,:,i], 1,r0,coords_evo[:,:,i],coords_evo_vertex[:,:,i])
    
    Force_center_vector[:,:,i] = F_center
    Force_vertex_vector[:,:,i] = F_vertex
    
    coords_evo[:,:,i+1] = coords_evo[:,:,i] + F_center*dt
    coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + F_vertex*dt
    vor_ridges_array[:,:,i+1], coords_evo_vertex[:,:,i+1], vor_regions =  sppv.T1transition(coords_evo_vertex[:,:,i+1],vor_ridges_array[:,:,i],vor_regions, vor_point_region,r0/2)
    
    vor_regions_array.append(vor_regions)
    
main_path = pathlib.Path().absolute()
datafilesavename = input('Name of savefile: ')
datafilesavename = datafilesavename + '.h5'

data_save = h5.File(datafilesavename,'w')
data_save.create_dataset('centers',data = coords_evo)
data_save.create_dataset('vertices',data=coords_evo_vertex)
data_save.create_dataset('forces_center',data=Force_center_vector)
data_save.create_dataset('forces_vertices',data=Force_vertex_vector)
data_save.create_dataset('Edge connections',data=vor_ridges_array)
data_save.create_dataset('number_of_points',data=N)
data_save.create_dataset('time_of_simul',data=M)
data_save.create_dataset('timestep',data = dt)
data_save.close()

