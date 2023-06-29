import numpy as np
import randomlatticethermalized as rlt
import selfpropelledparticlevoronoi as sppv
from scipy.spatial import Voronoi
import h5py as h5
import pathlib
import findconnections as fc

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
    N_m = int(input("Square root of number of points you want to add onto the network"))
    
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



    
vorRegions = []
vorPointRegion = []

for i in range(len(vorPointRegion1)):
    vorPointRegion.append(i)
    vorRegions.append(vorRegions1[vorPointRegion1[i]])

    
if (size_of_wound == 1) or (size_of_wound == 0):
    print('Single cell wound')
# else: 
#     wound_cells_loc = np.argsort(L_array)[:size_of_wound]
#     print(wound_cells_loc)
#     print(wound_loc)
#     #Removing the edge associated with the intersection between wound elements.
    
    
#     UnionIntersect = []
#     for i in range(size_of_wound):
#         for j in range(i+1,size_of_wound):
#             Ri = vorRegions[wound_cells_loc[i]]           
#             Rj = vorRegions[wound_cells_loc[j]]
#             UnionIntersect.append(list(set(Ri).intersection(Rj)))
#     print(UnionIntersect)   
#     for edge in UnionIntersect:
#         if len(edge)>1:
#             vorRidges.remove(sorted(edge))
        
    
#     #Removing the vertices associated with the intersection between 3 wound elements
#     if size_of_wound >= 3:
#         TripleIntersect = []
#         for i in range(size_of_wound):
#             for j in range(i+1,size_of_wound):
#                 for k in range(j+1, size_of_wound):
#                     Ri = vorRegions[wound_cells_loc[i]]          
#                     Rj = vorRegions[wound_cells_loc[j]]
#                     Rk = vorRegions[wound_cells_loc[k]]
                    
#                     TripleIntersect.append(list(set(Ri).intersection(Rj).intersection(Rk)))
#         print(TripleIntersect)
    
#         #vorVertex = np.delete(vorVertex,np.array(TripleIntersect),0)
                    
                    
        
#     #Merge the regions of neighbouring cells together and remove them from the list
    
#     wound_cells_loc= np.delete(wound_cells_loc,0)
    
#     print("region "+str(vorRegions[wound_loc])+", point = "+str(wound_loc))
#     for i in range(size_of_wound-1):
#         #for v in vorRegions[vorPointRegion[wound_cells_loc[i]]]:
#         for v in vorRegions[wound_cells_loc[i]]:
#             #if v not in vorRegions[vorPointRegion[wound_loc]]:
#             if v not in vorRegions[wound_loc]:
#                 #vorRegions[vorPointRegion[wound_loc]].append(v)
#                 vorRegions[wound_loc].append(v)
    
                
#     if size_of_wound >= 3:           
#         for vr in TripleIntersect:
#             if len(vr)> 0:
#                 #vorRegions[vorPointRegion[wound_loc]].remove(vr[0])
#                 vorRegions[wound_loc].remove(vr[0])
                
#     print("region "+str(vorRegions[wound_loc])+", point = "+str(wound_loc))
#     print(wound_cells_loc)
#     for i in range(size_of_wound-1):
#         vorRegions.pop(wound_cells_loc[i]) #remove elements from location

#     coords = np.delete(coords,np.array(wound_cells_loc),0)
    
    
N_new = len(coords)
L_array = np.zeros(N_new)
for i in range(N_new):
    L_array[i] = np.linalg.norm(coords[i])
    
wound_loc = np.argmin(L_array)
N = N_new
        
print("region "+str(vorRegions[wound_loc])+", point = "+str(wound_loc))

# Define arrays to insert numerical results

import matplotlib.pyplot as plt



#av = []

#for i in range(N):
#    av.append(sppv.area_vor(vorPointRegion,vorRegions,np.array(fc.vertices_in_bound(list(vor.vertices),5)),i))

#Maybe modify things a bit
L_max = 5
L_min = -5  
DeltaL = L_max - L_min



for i in range(N):
    with open('tissue225_1/tissue'+str(size_of_wound)+'/centers.txt','a') as text:
        text.write(str(i)+";"+str(coords[i,0])+";"+str(coords[i,1])+";"+str(vorPointRegion[i])+";"+str(vorRegions[vorPointRegion[i]])+"\n")



print("region "+str(vorRegions[wound_loc])+", point = "+str(wound_loc))

for v in range(len(vorVertex)):
    #print(v)
    NeighR, Neigh_C = fc.find_vertex_neighbour_centers(vorRegions,vorPointRegion,v)
    Neigh_V = fc.find_vertex_neighbour_vertices(vorRidges,v)
    with open('tissue225_1/tissue'+str(size_of_wound)+'/vertices.txt','a') as text:
        text.write(str(v)+";"+str(vorVertex[v,0])+";"+str(vorVertex[v,1])+";"+str(Neigh_C)+";"+str(Neigh_V)+"\n")
        
for e in range(len(vorRidges)):
    with open('tissue225_1/tissue'+str(size_of_wound)+'/edges.txt','a') as text:
        text.write(str(e)+";"+str(vorRidges[e])+"\n")

#First find the boundary
boundary_tissue = fc.find_boundary_vertices(len(vorVertex),vorRidges)
truebound = np.where(vorVertex[boundary_tissue,0]**2+vorVertex[boundary_tissue,1]**2>9)
boundary_tissue= list(np.array(boundary_tissue)[truebound])
boundary_wound = vorRegions[wound_loc]

with open('tissue225_1/tissue'+str(size_of_wound)+'/boundaries.txt','a') as text:
    text.write(str(boundary_tissue)+"\n")
    text.write(str(boundary_wound)+"\n")


with open('tissue225_1/tissue'+str(size_of_wound)+'/woundloc.txt','a') as text:
    text.write(str(wound_loc)+"\n")

print(wound_loc)
print("region "+str(boundary_wound))