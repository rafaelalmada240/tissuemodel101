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
print(np.mean(av))
print(np.std(av))
#Maybe modify things a bit
L_max = 5
L_min = -5  
DeltaL = L_max - L_min

r0 = np.sqrt(np.mean(av)/np.pi)


h = 1e-4*DeltaL/(2*np.sqrt(N)) #size step
print(h)

#for PC


#for cluster
#M = 50000 

#Model parameters
K_run = 1
A0_run = av#np.median(av)
G_run = 1
L_List = [0, 5, 10,15,20, 25, 30, 35, 40,45,50,55]
#[20, 21, 22, 23, 24, 25, 26, 27,28,29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,41, 42, 43, 44,45,46, 47,48,49,50,51,52,53,54,55]
# #float(input('Max Value of L/G ratio (times 10): '))
#Lw_run = float(input('Max Value of Lw ratio: '))
Lw_List = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

#[0.0, 0.2, 0.4, 0.6, 0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4, 2.8, 3.2, 3.6, 4.0]
mu = 1

dt = mu/(K_run*np.mean(av))*1e-2 #time step bigger than K_run A0_run^2/mu
print(dt)
T = 10
M = int(T/dt)
#coords_evo = np.zeros((N,2,M))
#coords_evo_vertex = np.zeros((len(vor.vertices),2,M))
#coords_evo[:,:,0] = coords
#coords_evo_vertex[:,:,0] = np.array(fc.vertices_in_bound(list(vor.vertices),2*L_max))





for Lr in L_List:#range(0,int(L_run)):
    # Run simulations
    #transition_counter = 0
    for Lw in Lw_List:
        
        coords_evo = coords
        coords_evo_vertex = np.array(fc.vertices_in_bound(list(vor.vertices),2*L_max))


        vorRidges = fc.remove_minus(vor.ridge_vertices)


        perimeterArray = np.zeros((coords_evo.shape[0],M-1))
        areaArray = np.zeros((coords_evo.shape[0],M-1))

        for i in range(M-1):
            #Compute forces
            F_vertex = sppv.force_vtx_elastic_wound(vorRegions, vorPointRegion, vorRidges, K_run,A0_run,G_run,-Lr/10,-Lw,coords_evo_vertex,coords_evo,wloc,h)
        
            #Compute the area and perimeter of the wound
    
            perimeterWound = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex,wloc)
            areaWound = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex,wloc)
        
            with open('simul3/woundinfoA'+str(np.median(A0_run).round(3))+'G'+str(G_run)+'L-'+str(Lr/10)+'Lw-'+str(Lw)+'Ncells'+str(N)+'.txt','a') as text:
                text.write(str(i))
                text.write(' ')
                text.write(str(perimeterWound))
                text.write(' ')
                text.write(str(areaWound)+ '\n')
                
            #Reflexive boundary conditions
            A_vertex = rlt.newWhere(coords_evo_vertex + mu*F_vertex*dt,6)
            coords_evo_vertex = coords_evo_vertex + mu*A_vertex*F_vertex*dt
        
            #Cell center positions are the average of the cell vertex positions
            coords_evo = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords_evo),np.array(coords_evo_vertex))
            vorRidges, coords_evo_vertex, vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex),np.array(vorRidges),vorRegions, vorPointRegion,0.01*r0)
        
            
