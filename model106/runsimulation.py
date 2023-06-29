import numpy as np
from scipy.spatial import Voronoi
import selfpropelledparticlevoronoi as sppv
import pathlib
import h5py as h5
import matplotlib.pyplot as plt
import randomlatticethermalized as rlt
import topologicaltransitions as tpt
import weightedgraphanalysis as wga
import findconnections as fc

''' 
Run simulations in this module for a prebuilt tissue network 
(that is we generate the tissue centers in a separate file, and we generate the vertices through a voronoi tesselation in this file)
'''

main_path = pathlib.Path().absolute()
datafileloadname = input('Name of savefile to load: ')
datafileloadname = datafileloadname + '.h5'
data_set = h5.File(str(main_path) + '/' + datafileloadname,'r')
coords = np.array(data_set['centers'])
N = np.array(data_set['number of centers'])
wloc = np.array(data_set['woundlocation'])
data_set.close()

# Define arrays to insert numerical results

vor = Voronoi(coords)
vorPointRegion = vor.point_region
vorRegions = vor.regions

av = []

for i in range(N):
    av.append(sppv.area_vor(vorPointRegion,vorRegions,np.array(fc.vertices_in_bound(list(vor.vertices),5)),i))
print(np.median(av))
#Maybe modify things a bit
L_max = 5
L_min = -5  
DeltaL = L_max - L_min

r0 = np.sqrt(np.mean(av)/np.pi)

dt = 2e-2 #time step
h = 1e-4*DeltaL/(2*np.sqrt(N)) #size step

#for PC
M = 3

#for cluster
#M = 50000 

#Model parameters
K_run = 1
A0_run = av#np.median(av)
G_run = 1
L_run = -float(input('Value of L/G ratio: '))
Lw_run = -float(input('Value of Lw ratio: '))
mu = 1


coords_evo = np.zeros((N,2,M))
coords_evo_vertex = np.zeros((len(vor.vertices),2,M))
coords_evo[:,:,0] = coords
coords_evo_vertex[:,:,0] = np.array(fc.vertices_in_bound(list(vor.vertices),2*L_max))

magF = float(input('Magnitude of force acting on the tissue: '))


vorRidges = fc.remove_minus(vor.ridge_vertices)
Force_vertex_vector = np.zeros((len(vor.vertices),2,M))

Adjacency_Mat = np.zeros((len(vor.vertices),len(vor.vertices),M))

vorRidges_array = np.zeros((len(vorRidges),2,M))
vorRidges_array[:,:,0] = vorRidges

vorRegions_array = []
vorRegions_array.append(vorRegions)

perimeterArray = np.zeros((coords_evo.shape[0],M-1))
areaArray = np.zeros((coords_evo.shape[0],M-1))

continue_option = int(input('Continue with simulation (y-1/n-0): '))
wound_option = int(input('Are you simulating a wound? (y-1/n-0): '))


if continue_option == 1:
    # Run simulations
    transition_counter = 0
    
    if wound_option == 1:

        for i in range(M-1):
            #Compute forces
            F_vertex = sppv.force_vtx_elastic_wound(vorRegions, vorPointRegion, vorRidges_array[:,:,i], K_run,A0_run,G_run,L_run,Lw_run,coords_evo_vertex[:,:,i],coords_evo[:,:,i],wloc,h)
            Force_vertex_vector[:,:,i] = F_vertex
        
            #Compute the area and perimeter for all cell elements in the network
            for elem in range(len(coords_evo)):
                perimeterArray[elem,i] = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
                areaArray[elem,i] = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
        
            with open('woundareaG'+str(G_run)+'L'+str(L_run)+'Lw'+str(Lw_run)+'Ncells'+str(N)+'.txt','a') as text:
                text.write(str(i))
                text.write(' ')
                text.write(str(areaArray[wloc,i])+ '\n')
            #Reflexive boundary conditions
            A_vertex = rlt.newWhere(coords_evo_vertex[:,:,i] + mu*F_vertex*dt,6)
            coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + mu*A_vertex*F_vertex*dt
        
            #Cell center positions are the average of the cell vertex positions
            coords_evo[:,:,i+1] = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords_evo[:,:,i]),np.array(coords_evo_vertex[:,:,i+1]))
        
            vorRidges_array[:,:,i+1] = vorRidges_array[:,:,i]
        
            #Do topological transition T1 (hopefully works)
            #vorRidges_array[:,:,i+1], coords_evo_vertex[:,:,i+1], vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex[:,:,i+1]),np.array(vorRidges_array[:,:,i]),vorRegions, vorPointRegion,0.01*r0)
        
            #Compute the adjacency matrix of the network
            Adjacency_Mat[:,:,i] = wga.adjacency_matrix(len(vor.vertices),np.array(vorRidges_array[:,:,i+1]))
        
            if transition_counter > 0:
                with open('transition_times.txt','a') as text:
                    text.write(str(i))
                    text.write(' ')
                    text.write(str(transition_counter)+ '\n')
        
            vorRegions_array.append(vorRegions)
            
    else:
        shear_option = int(input('Are you simulating a shear of the tissue? (y-1/n-0): '))
        if shear_option ==1: 
            
            #smooth boundaries
            
            vorVertices = np.array(fc.vertices_in_bound(list(vor.vertices),L_max))
            bound = fc.find_boundary_vertices(len(vorVertices),vorRidges)
            
            x_down = -3
            x_up = 3
            lbound_down = np.where(vorVertices[bound,1]<x_down)[0]
            bound_down = np.array(bound)[lbound_down]
            lbound_up = np.where(vorVertices[bound,1]>x_up)[0]
            bound_up = np.array(bound)[lbound_up]
            #
            
            for i in range(30):
                #Compute forces
                F_vertex = sppv.smooth_boundary(vorRegions, vorPointRegion, vorRidges, K_run,A0_run,G_run,L_run,vorVertices,coords,h)
                
                
                #Reflexive boundary conditions
                A_vertex = rlt.newWhere(coords_evo_vertex[:,:,i] + mu*F_vertex*dt,10)
                vorVertices = vorVertices + mu*A_vertex*F_vertex*dt
        
                #Cell center positions are the average of the cell vertex positions
                coords = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords),np.array(vorVertices))
        
            #Initial relaxation
            
            N_relax = int(M/10)
            for i in range(N_relax-1):
                #Compute forces
                F_vertex = sppv.force_vtx_elastic_shear_rest(vorRegions, vorPointRegion, vorRidges, K_run,A0_run,G_run,L_run,vorVertices,coords,h,bound, bound_down)
                
                #Reflexive boundary conditions
                A_vertex = rlt.newWhere(coords_evo_vertex[:,:,i] + mu*F_vertex*dt,6)
                vorVertices = vorVertices + mu*A_vertex*F_vertex*dt
        
                #Cell center positions are the average of the cell vertex positions
                coords = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords),np.array(vorVertices))
        
                #vorRidges = vorRidges
        
                #Do topological transition T1 (hopefully works)
                vorRidges, vorVertices, vorRegions, transition_counter =  tpt.T1transition2(np.array(vorVertices),np.array(vorRidges),vorRegions, vorPointRegion,0.01*r0)
              
            coords_evo[:,:,0] = coords
            coords_evo_vertex[:,:,0] = np.array(fc.vertices_in_bound(list(vorVertices),L_max))
            vorRidges_array[:,:,0] = vorRidges
            vorRegions_array = []
            vorRegions_array.append(vorRegions)
        
            for i in range(M-1):
                #Compute forces
                F_vertex = sppv.force_vtx_elastic_shear(vorRegions, vorPointRegion, vorRidges_array[:,:,i], K_run,A0_run,G_run,L_run,coords_evo_vertex[:,:,i],coords_evo[:,:,i],magF,h, bound, bound_up, bound_down)
                Force_vertex_vector[:,:,i] = F_vertex
        
                #Compute the area and perimeter for all cell elements in the network
                for elem in range(len(coords_evo)):
                    perimeterArray[elem,i] = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
                    areaArray[elem,i] = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
                
                #Reflexive boundary conditions
                #A_vertex = rlt.newWhere(coords_evo_vertex[:,:,i] + mu*F_vertex*dt,6)
                coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + mu*F_vertex*dt
        
                #Cell center positions are the average of the cell vertex positions
                coords_evo[:,:,i+1] = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords_evo[:,:,i]),np.array(coords_evo_vertex[:,:,i+1]))
        
                #vorRidges_array[:,:,i+1] = vorRidges_array[:,:,i]
        
                #Do topological transition T1 (hopefully works)
                vorRidges_array[:,:,i+1], coords_evo_vertex[:,:,i+1], vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex[:,:,i+1]),np.array(vorRidges_array[:,:,i]),vorRegions, vorPointRegion,0.01*r0)
        
                #Compute the adjacency matrix of the network
                Adjacency_Mat[:,:,i] = wga.adjacency_matrix(len(vor.vertices),np.array(vorRidges_array[:,:,i+1]))
        
                if transition_counter > 0:
                    with open('transition_times.txt','a') as text:
                        text.write(str(i))
                        text.write(' ')
                        text.write(str(transition_counter)+ '\n')
        
                vorRegions_array.append(vorRegions)
            
        else:
            for i in range(M-1):
                #Compute forces
                F_vertex = sppv.force_vtx_elastic(vorRegions, vorPointRegion, vorRidges_array[:,:,i], K_run,A0_run,G_run,L_run,coords_evo_vertex[:,:,i],coords_evo[:,:,i],h)
                Force_vertex_vector[:,:,i] = F_vertex
        
                #Compute the area and perimeter for all cell elements in the network
                for elem in range(len(coords_evo)):
                    perimeterArray[elem,i] = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
                    areaArray[elem,i] = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex[:,:,i],elem)
                
                #Reflexive boundary conditions
                A_vertex = rlt.newWhere(coords_evo_vertex[:,:,i] + mu*F_vertex*dt,6)
                coords_evo_vertex[:,:,i+1] = coords_evo_vertex[:,:,i] + mu*A_vertex*F_vertex*dt
        
                #Cell center positions are the average of the cell vertex positions
                coords_evo[:,:,i+1] = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords_evo[:,:,i]),np.array(coords_evo_vertex[:,:,i+1]))
        
                #vorRidges_array[:,:,i+1] = vorRidges_array[:,:,i]
        
                #Do topological transition T1 (hopefully works)
                vorRidges_array[:,:,i+1], coords_evo_vertex[:,:,i+1], vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex[:,:,i+1]),np.array(vorRidges_array[:,:,i]),vorRegions, vorPointRegion,0.3*r0)
        
                #Compute the adjacency matrix of the network
                Adjacency_Mat[:,:,i] = wga.adjacency_matrix(len(vor.vertices),np.array(vorRidges_array[:,:,i+1]))
        
                if transition_counter > 0:
                    with open('transition_times.txt','a') as text:
                        text.write(str(i))
                        text.write(' ')
                        text.write(str(transition_counter)+ '\n')
        
                vorRegions_array.append(vorRegions)
    

    
   # Save data in a file
            
    main_path = pathlib.Path().absolute()
    datafilesavename = input('Name of savefile: ')
    datafilesavename = datafilesavename + '.h5'
    data_save = h5.File(datafilesavename,'w')
    data_save.create_dataset('centers',data = coords_evo)
    data_save.create_dataset('vertices',data=coords_evo_vertex)
    data_save.create_dataset('forces_vertices',data=Force_vertex_vector)
    data_save.create_dataset('Edge connections',data=vorRidges_array)
    data_save.create_dataset('AdjacencyMatrix',data=Adjacency_Mat)
    data_save.create_dataset('perimeter',data=perimeterArray)
    data_save.create_dataset('area',data=areaArray)
    data_save.create_dataset('wound area', data = areaArray[wloc])
    data_save.create_dataset('number_of_points',data=N)
    data_save.create_dataset('time_of_simul',data=M)
    data_save.create_dataset('timestep',data = dt)
    data_save.close()