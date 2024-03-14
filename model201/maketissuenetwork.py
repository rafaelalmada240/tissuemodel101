import numpy as np
import randomlatticethermalized as rlt
import selfpropelledparticlevoronoi as sppv
from scipy.spatial import Voronoi
import findconnections as fc
import itertools
import readTissueFiles as rTF
import time

ti = time.time()
def newAbs(a,b):
    res = np.zeros(a.shape)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if np.abs(a[i,j])<b:
                res[i,j] = a[i,j]
            else:
                res[i,j] = np.sign(a[i,j])*b
    return res 

option_1 = int(input(" Do you want to generate a random lattice: y(1) or n(0): "))

if option_1 == 0:
    # Define a regular lattice
    option_2 = int(input("Use which method to generated regular lattice: "))
    N_m = int(input("Square root of number of points you want to add onto the network"))
    
    if option_2 == 0:
    
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
    else:
        # Generate a regular lattice in a square domain
        L_max = 5
        L_min = -5  
        dL = (L_max-L_min)/N_m
        A_0 = np.array([[1,np.sqrt(1)/2],[0,1]])
        u = np.linspace(L_min,L_max,N_m)
        v = np.linspace(L_min,L_max,N_m)

        coord_list = []

        for j in range(N_m):
            for i in range(N_m):
                coord_list.append(A_0.dot(np.array([u[j],v[i]])))
            
        coords = np.array(coord_list)
        N = N_m**2
        
        for i in range(N):
            if coords[i,0]<-5:
                coords[i,0] = coords[i,0]+dL*N_m+dL+np.sqrt(1)/(2*N_m)
            if coords[i,0] > 5:
                coords[i,0] = coords[i,0]-dL*N_m-dL-np.sqrt(1)/(2*N_m)
            


if option_1 == 1:
    # Define a random lattice

    N = int(input("Square root of the number of points you want to add onto the network "))**2
    coords = newAbs(rlt.newcoords(N),5)
    

  
# For wound modelling, find a cell region to classify as a wound
#An initial approach is to choose the region closest to the center

L_array = np.zeros(N)
for i in range(N):
    L_array[i] = np.linalg.norm(coords[i])
    
wound_loc = np.argmin(L_array)
L_max = 5
L_min = -5  
dL = (L_max-L_min)/N**0.5
size_of_wound = int(input("Number of cells to replace with wound: "))

vor = Voronoi(coords)
vorPointRegion1 = vor.point_region
vorRegions1 = vor.regions
vorVertex = np.array(fc.vertices_in_bound(list(vor.vertices),2*L_max))
vorRidges = fc.remove_minus(vor.ridge_vertices)


for e in range(len(vorRidges)):
    vorRidges[e] = sorted(vorRidges[e])
    
vorRegions = []
vorPointRegion2 = []

for i in range(len(vorPointRegion1)):
    vorPointRegion2.append(i)
    vorRegions.append(vorRegions1[vorPointRegion1[i]])
    
if (size_of_wound == 1) or (size_of_wound == 0):
    print('Single cell wound')

else: 
    wound_cells_loc = np.argsort(L_array)[:size_of_wound]
    print('Multiple cell wounds at: ')

    #Removing the edge associated with the intersection between wound elements.
    
    UnionIntersect = []
    List_pair = itertools.permutations(wound_cells_loc,2)
    for listit in List_pair:
        Ri = vorRegions[listit[0]]           
        Rj = vorRegions[listit[1]]
        eij = sorted(list(set(Ri).intersection(Rj)))
        if (eij not in UnionIntersect) and (len(eij)>0):
            UnionIntersect.append(eij)
    print("Edge shared by the regions, to be removed: ")

    for edge in UnionIntersect:
        if len(edge)>1:
            vorRidges.remove(sorted(edge))

    #Removing the vertices associated with the intersection between 3 wound elements
    if size_of_wound >= 3:
        List_iter = itertools.permutations(wound_cells_loc, 3)
        TripleIntersect = []
        for listiter in List_iter:
            Ri = vorRegions[listiter[0]]          
            Rj = vorRegions[listiter[1]]
            Rk = vorRegions[listiter[2]]
            vijk = list(set(Ri).intersection(Rj).intersection(Rk))
            if (vijk not in TripleIntersect) and (len(vijk)>0):
                TripleIntersect.append(vijk)
        flat_list = []
        for row in TripleIntersect:
            flat_list += row
        TripleIntersect = flat_list
        print("Vertex shared by multiple regions, to be removed: ")
        
    #Merge the regions of neighbouring cells together and remove them from the list
    
    wound_cells_loc= np.delete(wound_cells_loc,0)
    print("Merging the following regions")


    for i in range(size_of_wound-1):
        for v in vorRegions[wound_cells_loc[i]]:
            if v not in vorRegions[wound_loc]:
                vorRegions[wound_loc].append(v)
    
                
    if size_of_wound >= 3:           
        for vr in TripleIntersect:
            for r in range(len(vorRegions)):
                if vr in vorRegions[r]:
                    vorRegions[r].remove(vr)
        # Removing isolated vertices in an edge
        
        
        # for v in vorRegions[wound_loc]:
        #     neighv = fc.find_vertex_neighbour_vertices(vorRidges,v)
        #     if len(neighv)<3:
        #         for vi in neighv:
                
                
        # The removal of vertices throws some labels out of whack. Maybe this will solve some of the issues 
    # TpInter = list(np.sort(TripleIntersect)[::-1])
    # print(TpInter)
    # nvorRegions = list(vorRegions)
    # nvorRidges = list(vorRidges)
    # for v in TpInter:
    #     print(v)
    #     for r in range(len(vorRegions)):
    #         for vr in range(len(vorRegions[r])):
    #             if (nvorRegions[r][vr] >= v) and (nvorRegions[r][vr] not in TpInter) :
    #                 vorRegions[r][vr] = vorRegions[r][vr] -1
    #     for e in range(len(vorRidges)):
    #         for ve in range(len(vorRidges[e])):
    #             if (nvorRidges[e][ve] >= v) and (nvorRidges[e][ve] not in TpInter):
    #                 vorRidges[e][ve] = vorRidges[e][ve] -1
                

    print("Removing intersect vertices")       

    nvorRegions = list(vorRegions)
    for i in range(size_of_wound-1):
        vorRegions.remove(nvorRegions[wound_cells_loc[i]]) #remove elements from location
    coords = np.delete(coords,np.array(wound_cells_loc),0)
    

    
# print(N)
N_new = len(coords)
L_array = np.zeros(N_new)
for i in range(N_new):
    L_array[i] = np.linalg.norm(coords[i])
    
wound_loc = np.argmin(L_array)
N = N_new

# Define arrays to insert numerical results

vorPointRegion = []

for i in range(N):
    vorPointRegion.append(i)
av = sppv.areas_vor(vorPointRegion,vorRegions,vorVertex,vorRidges,vorPointRegion)

#Maybe modify things a bit
L_max = 5
L_min = -5  
DeltaL = L_max - L_min

foldername='tissue'
rTF.saveGeneratedNetwork(coords, size_of_wound,vorPointRegion,vorRegions,vorVertex,vorRidges,wound_loc,foldername)

tf = time.time() - ti
with open(foldername+'/log_tissue_generation'+str(size_of_wound)+'.txt','a') as log:
    
    log.write('Random Lattice? - ')
    log.write(' ')
    log.write(str(option_1)+'\n')
    
    if size_of_wound== 1:
        log.write('Single cell wound at - ')
        log.write(' ')
        log.write(str(wound_loc)+'\n')
    else:    
        log.write('Multiple cells wound at - ')
        log.write(' ')
        log.write(str(wound_loc)+'\n')
        
    log.write('Edge shared by the regions, to be removed - ')
    log.write(' ')
    log.write(str(UnionIntersect)+'\n')

    log.write('Vertex shared by multiple regions, to be removed - ')
    log.write(' ')
    log.write(str(TripleIntersect)+'\n')
    
    log.write('Regions removed - ')
    log.write(' ')
    log.write(str(wound_cells_loc)+'\n')
    
    log.write('Number of remaining regions - ')
    log.write(' ')
    log.write(str(N_new)+'\n')
    
    log.write('Simulation time (s) - ')
    log.write(' ')
    log.write(str(tf.round(3))+'\n')
