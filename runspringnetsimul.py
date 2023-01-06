import numpy as np
from scipy.spatial import Voronoi
import selfpropelledparticlevoronoi as sppv
import pathlib
import h5py as h5

''' 

Run simulations in this module

'''
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

Force_center_vector = np.zeros((N,2,M))
Force_vertex_vector = np.zeros((len(vor.vertices),2,M))

# Run simulations

for i in range(M-1):
 
    F_center, F_vertex = sppv.force_vor_elastic(vor,1,r0,coords_evo[:,:,i],coords_evo_vertex[:,:,i])
    
    Force_center_vector[:,:,i] = F_center
    Force_vertex_vector[:,:,i] = F_vertex
    
    coords_evo[:,:,i+1] = coords_evo[:,:,i] + F_center*dt
    coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + F_vertex*dt
    
main_path = pathlib.Path().absolute()
datafilesavename = input('Name of savefile: ')
datafilesavename = datafilesavename + '.h5'

data_save = h5.File(datafilesavename,'w')
data_save.create_dataset('centers',data = coords_evo)
data_save.create_dataset('vertices',data=coords_evo_vertex)
data_save.create_dataset('forces_center',data=Force_center_vector)
data_save.create_dataset('forces_vertices',data=Force_vertex_vector)
data_save.create_dataset('Edge connections',data=vor.ridge_vertices)
data_save.create_dataset('number_of_points',data=N)
data_save.create_dataset('time_of_simul',data=M)
data_save.create_dataset('timestep',data = dt)
data_save.close()

