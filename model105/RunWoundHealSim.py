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

wound_sizes = [2,3,4,5,6]

for size_of_wound in wound_sizes:
    
    print("wound size: "+str(size_of_wound))

    coords = []
    vorPointRegion1 = []
    vorRegions = []

    with open('closurewithtransition225/size'+str(size_of_wound)+'/centers.txt','r') as text:
                for line in text:
                    last_line = (line.replace("\n","")).split(';')
                    coords.append([float(last_line[1]),float(last_line[2])])
                    vorPointRegion1.append(int(last_line[3]))
                    l4 = last_line[4].replace("[","").replace("]","").split(',')
                    lint4 = []
                    if all('' == s or s.isspace() for s in l4):
                        vorRegions.append(lint4)
                        continue
                    else:
                        for r in l4:
                            lint4.append(int(r))
                    vorRegions.append(lint4)
                    
    coords = np.array(coords)
    N = len(coords[:,0])
    print(N)

    vertices = []
    with open('closurewithtransition225/size'+str(size_of_wound)+'/vertices.txt','r') as text:
                for line in text:
                    last_line = (line.replace("\n","")).split(';')
                    vertices.append([float(last_line[1]),float(last_line[2])])

    #print(vorPointRegion)

    vorRidges = []
    with open('closurewithtransition225/size'+str(size_of_wound)+'/edges.txt','r') as text:
                for line in text:
                    last_line = (line.replace("\n","")).split(';')
                    l4 = last_line[1].replace("[","").replace("]","").split(',')
                    lint4 = []
                    if all('' == s or s.isspace() for s in l4):
                        vorRidges.append(lint4)
                        continue
                    else:
                        for r in l4:
                            lint4.append(int(r))
                    vorRidges.append(lint4)
    Boundaries = []                
    with open('closurewithtransition225/size'+str(size_of_wound)+'/boundaries.txt','r') as text:
        for line in text:
            last_line = line.replace("\n","")
            l4 = line.replace("[","").replace("]","").split(',')
            lint4 = []
            if all('' == s or s.isspace() for s in l4):
                Boundaries.append(lint4)
                continue
            else:
                for r in l4:
                    lint4.append(int(r))
            Boundaries.append(lint4)

    wloc = 0

    with open('closurewithtransition225/size'+str(size_of_wound)+'/woundloc.txt','r') as text:
        for line in text:
            wloc = int(line.replace("\n",""))
        
    av = []
    
    vertices = np.array(vertices)

    vorPointRegion= []
    for i in range(N):
        vorPointRegion.append(i)

    for i in range(N):
        av.append(sppv.area_vor(vorPointRegion,vorRegions,vertices,vorRidges,i))
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
    L_List =  [0, .5, 1.0,1.5,2.0, 2.5, 3.0, 3.5, 4.0,4.5,5.0,5.5]
    #[20, 21, 22, 23, 24, 25, 26, 27,28,29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,41, 42, 43, 44,45,46, 47,48,49,50,51,52,53,54,55]
    # #float(input('Max Value of L/G ratio (times 10): '))
    #Lw_run = float(input('Max Value of Lw ratio: '))
    Lw_List = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    #[0.0, 0.2, 0.4, 0.6, 0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4, 2.8, 3.2, 3.6, 4.0]
    mu = 1

    dt = mu/(K_run*np.mean(av))*25e-4 #time step bigger than K_run A0_run^2/mu
    print(dt)
    T = 10
    M = int(T/dt)
    
    areaWound0 = sppv.area_vor(vorPointRegion,vorRegions,vertices,vorRidges,wloc)
    
    ordered_boundary = [0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 0]
    
    


    
    for Lr in L_List:#range(0,int(L_run)):
        # Run simulations
        #transition_counter = 0
        for Lw in Lw_List:
            
            i = 0
            areaWound = areaWound0
            coords_evo = coords
            coords_evo_vertex = vertices

            while (i < M) and ((areaWound>= areaWound0/2) and (areaWound<= 2*areaWound0)):
                #Compute forces
                F_vertex = sppv.force_vtx_elastic_wound(vorRegions, vorPointRegion, vorRidges, K_run,A0_run,G_run,-Lr,-Lw,coords_evo_vertex,coords_evo,wloc,h)
                #print(i)
            
                #Compute the area and perimeter of the wound
        
                perimeterWound = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex,vorRidges,wloc)
                areaWound = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex,vorRidges,wloc)
            
                
                    
                #Reflexive boundary conditions
                A_vertex = rlt.newWhere(coords_evo_vertex + mu*F_vertex*dt,6)
                coords_evo_vertex = coords_evo_vertex + mu*A_vertex*F_vertex*dt
                
                
            
                #Cell center positions are the average of the cell vertex positions
                coords_evo = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords_evo),np.array(coords_evo_vertex))
                
                #vorRidges, coords_evo_vertex, vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex),np.array(vorRidges),vorRegions, vorPointRegion,0.01*r0)
                with open('closurewithtransition225/size'+str(size_of_wound)+'/woundinfoA'+str(areaWound0.round(3))+'G'+str(G_run)+'L-'+str(Lr)+'Lw-'+str(Lw)+'Ncells'+str(N)+'.txt','a') as text:
                    text.write(str(i))
                    #text.write(' ')
                    #text.write(str(transition_counter))
                    text.write(' ')
                    text.write(str(perimeterWound))
                    text.write(' ')
                    text.write(str(areaWound)+ '\n')
                    
                i = i + 1
                
                
                f_vert = fc.norm(F_vertex)
                n_vert = f_vert/np.max(f_vert)
                
                #plt.figure(figsize=(6,6))
                #plt.plot(coords_evo_vertex[list(np.array(Boundaries[1])[ordered_boundary]),0],coords_evo_vertex[list(np.array(Boundaries[1])[ordered_boundary]),1],'ro-')
                #plt.plot(coords_evo_vertex[:,0],coords_evo_vertex[:,1],'o', alpha = n_vert)
                #plt.grid()
                #plt.xlim(-6,6)
                #plt.ylim(-6,6)
                #plt.savefig('closurewithtransition225/size'+str(size_of_wound)+'/AverageDisplacementsForces'+str(i)+'.png',dpi=150)
                #plt.close()
            
                
