import numpy as np
from scipy.spatial import Voronoi
import selfpropelledparticlevoronoi as sppv
import pathlib
import h5py as h5
import matplotlib.pyplot as plt
import randomlatticethermalized as rlt
import topologicaltransitions as tpt
''' 

Run simulations in this module

'''

def newAbs(a,b):
    res = np.zeros(a.shape)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if np.abs(a[i,j])<b:
                res[i,j] = a[i,j]
            else:
                res[i,j] = np.sign(a[i,j])*b
    return res 

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

    coords = newAbs(rlt.newcoords(N),5)


# Define arrays to insert numerical results

L_max = 5
L_min = -5  
dt = 1e-3

M = 100

vor = Voronoi(coords)

vorPointRegion = vor.point_region
vorRegions = vor.regions

av = []

for i in range(N):
    av.append(sppv.area_vor(vorPointRegion,vorRegions,vor.vertices,i))

DeltaL = L_max - L_min
#r0 = DeltaL/(0.66*np.sqrt(N))
r0 = np.sqrt(np.median(av)/np.pi)
print(DeltaL/(r0))
#Model parameters
K_run = 0.001#0.01#0*0.05
A0_run = np.pi*r0**2
#P0 = 3.0
G_run = 0.15#10*K_run#0*0.5+0*0.25*A0_run*K_run
L_run = -0.25#G_run*P0*K_run*np.sqrt(A0_run)#-0*1-0*1*K_run*A0_run**(3/2)


coords_evo = np.zeros((N,2,M))
coords_evo_vertex = np.zeros((len(vor.vertices),2,M))
coords_evo[:,:,0] = coords
coords_evo_vertex[:,:,0] = newAbs(vor.vertices,5)

#coords_evo_v = np.zeros((N,2,M))
#coords_evo_vertex_v = np.zeros((len(vor.vertices),2,M))




vorRidges = sppv.remove_minus(vor.ridge_vertices)

#Force_center_vector = np.zeros((N,2,M))
Force_vertex_vector = np.zeros((len(vor.vertices),2,M))

vorRidges_array = np.zeros((len(vorRidges),2,M))
vorRidges_array[:,:,0] = vorRidges

vorRegions_array = []
vorRegions_array.append(vorRegions)

perimeterArray = np.zeros((coords_evo.shape[0],M-1))
areaArray = np.zeros((coords_evo.shape[0],M-1))
#boundary = sppv.find_boundary_vertices(len(vor.vertices),vorRidges)

continue_option = int(input('Continue with simulation: (y-1/n-0): '))

if continue_option == 1:
    
    # Run simulations
    transition_counter = 0

    for i in range(M-1):
        #Do T1 transitions and compute forces afterwards

        F_vertex = sppv.force_vtx_elastic(vorRegions, vorPointRegion, vorRidges_array[:,:,i], K_run,A0_run,G_run,L_run,coords_evo_vertex[:,:,i],coords_evo[:,:,i],0.01)
        Force_vertex_vector[:,:,i] = F_vertex
        
        for elem in range(len(coords_evo)):
            perimeterArray[elem,i] = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
            areaArray[elem,i] = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
        #Reflexive boundary conditions
        
        A_vertex = rlt.newWhere(coords_evo_vertex[:,:,i] + F_vertex*dt,6)
        coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + A_vertex*F_vertex*dt
        
        #Cell center positions are the average of the cell vertex positions
        coords_evo[:,:,i+1] = sppv.cells_avg_vtx(vorRegions,vorPointRegion,coords_evo[:,:,i],coords_evo_vertex[:,:,i+1])
        
        vorRidges_array[:,:,i+1] = vorRidges_array[:,:,i]
        
        #Do topological transition T1 (hopefully works)
        #vorRidges_array[:,:,i+1], coords_evo_vertex[:,:,i+1], vorRegions, transition_counter =  tpt.T1transition(coords_evo_vertex[:,:,i+1],vorRidges_array[:,:,i],vorRegions, vorPointRegion,r0/4)
        
        if transition_counter > 0:
            with open('transition_times.txt','a') as text:
                text.write('At time ')
                text.write(str(i))
                text.write(' number of transitions: ')
                text.write(str(transition_counter)+ '\n')
        
        vorRegions_array.append(vorRegions)
    
    #vorRegions_array = np.array(vorRegions_array)
    
   
            
    main_path = pathlib.Path().absolute()
    datafilesavename = input('Name of savefile: ')
    datafilesavename = datafilesavename + '.h5'

    data_save = h5.File(datafilesavename,'w')
    data_save.create_dataset('centers',data = coords_evo)
    data_save.create_dataset('vertices',data=coords_evo_vertex)
    #data_save.create_dataset('forces_center',data=Force_center_vector)
    data_save.create_dataset('forces_vertices',data=Force_vertex_vector)
    data_save.create_dataset('Edge connections',data=vorRidges_array)
    #data_save.create_dataset('regions',data=vorRegions_array, dtype='f4')
    data_save.create_dataset('perimeter',data=perimeterArray)
    data_save.create_dataset('area',data=areaArray)
    #data_save.create_dataset('point_regions',data=np.array(vorPointRegion))
    data_save.create_dataset('number_of_points',data=N)
    data_save.create_dataset('time_of_simul',data=M)
    data_save.create_dataset('timestep',data = dt)
    data_save.close()

