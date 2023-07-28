import numpy as np
import selfpropelledparticlevoronoi as sppv

import randomlatticethermalized as rlt
import topologicaltransitions as tpt
import findconnections as fc

''' 
Run simulations in this module for a prebuilt tissue network 
'''

size_of_wound = 1

#Loading up the dataset for the generated tissue

print("wound size: "+str(size_of_wound))

coords = []
vorPointRegion1 = []
vorRegions = []

with open('tissue/size'+str(size_of_wound)+'/centers.txt','r') as text:
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
with open('tissue/size'+str(size_of_wound)+'/vertices.txt','r') as text:
            for line in text:
                last_line = (line.replace("\n","")).split(';')
                vertices.append([float(last_line[1]),float(last_line[2])])


vorRidges = []
with open('tissue/size'+str(size_of_wound)+'/edges.txt','r') as text:
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
with open('tissue/size'+str(size_of_wound)+'/boundaries.txt','r') as text:
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

with open('tissue/size'+str(size_of_wound)+'/woundloc.txt','r') as text:
    for line in text:
        wloc = int(line.replace("\n",""))
    
av = []



vertices = np.array(vertices)

vorPointRegion= []
for i in range(N):
    vorPointRegion.append(i)


av = sppv.areas_vor(vorPointRegion,vorRegions,vertices,vorRidges,vorPointRegion)
#for i in range(N):
#    av.append(sppv.area_vor(vorPointRegion,vorRegions,vertices,vorRidges,i))
print(np.mean(av))
print(np.std(av))


#Unfortunate naming (L here means length)
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

#Instead of going with a full list, you may try initially with some parameters
#(L here stands for Lambda and represents linear tension in the energy functional)
Lr = float(input('Value of L/G ratio: '))
Lw = float(input(' Value of Lw ratio: '))

mu = 1

dt = mu/(K_run*np.mean(av))*25e-4 #time step bigger than K_run A0_run^2/mu
print(dt)
T = 10
M = int(T/dt)

areaWound0 = sppv.area_vor(vorPointRegion,vorRegions,vertices,vorRidges,wloc)

simple_output = int(input('Do you want a simple output? (y-1/n-0): '))

f_imp = 0
sigma_r = 0.0
i = 0
areaWound = areaWound0
coords_evo = coords
coords_evo_vertex = vertices

while (i < M) and ((areaWound>= areaWound0/4) and (areaWound<= 4*areaWound0)):
    #Compute forces
    F_vertex = sppv.force_vtx_elastic_wound(vorRegions, vorPointRegion, vorRidges, K_run,A0_run,G_run,Lr,Lw,coords_evo_vertex,coords_evo,wloc,f_imp,h)
    
    boundary_wound = fc.find_wound_boundary(vorRegions,vorPointRegion,wloc)
    
    
    
    
    for vf in boundary_wound:
    	#If direction of impulse is random for each vertex, use value of sigma diff from 0
    	#Random directional diffusion with drift
    	v_imp = np.array([1,0])+sigma_r*(np.random.rand(2)*2-1)
    	n_imp = v_imp/fc.norm(v_imp)
    	F_vertex[vf] = F_vertex[vf] + f_imp*n_imp
    #print(i)

    #Compute the area and perimeter of the wound

    perimeterWound = sppv.perimeter_vor(vorPointRegion,vorRegions,coords_evo_vertex,vorRidges,wloc)
    areaWound = sppv.area_vor(vorPointRegion,vorRegions,coords_evo_vertex,vorRidges,wloc)

    
        
    #Reflexive boundary conditions
    A_vertex = rlt.newWhere(coords_evo_vertex + mu*F_vertex*dt,6)
    coords_evo_vertex = coords_evo_vertex + mu*A_vertex*F_vertex*dt
    
    

    #Cell center positions are the average of the cell vertex positions
    coords_evo = sppv.cells_avg_vtx(vorRegions,vorPointRegion,np.array(coords_evo),np.array(coords_evo_vertex))
    
    vorRidges, coords_evo_vertex, vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex),np.array(vorRidges),vorRegions, vorPointRegion,0.01*r0)
    if simple_output == 1:
        with open('simpleoutput/woundinfoA'+str(areaWound0.round(3))+'G'+str(G_run)+'L-'+str(Lr)+'Lw-'+str(Lw)+'Ncells'+str(N)+'.txt','a') as text:
            text.write(str(i))
            text.write(' ')
            text.write(str(transition_counter))
            text.write(' ')
            text.write(str(coords_evo[wloc,0]))
            text.write(' ')
            text.write(str(coords_evo[wloc,1]))
            text.write(' ')
            text.write(str(perimeterWound))
            text.write(' ')

            text.write(str(areaWound)+ '\n')
    if simple_output == 0:
        #Need to change this
        for c in range(len(coords_evo)):
            with open('movieoutput/centers'+str(i)+'.txt','a') as text:
                text.write(str(c)+";"+str(coords_evo[c,0])+";"+str(coords_evo[c,1])+";"+str(vorPointRegion[c])+";"+str(vorRegions[vorPointRegion[c]])+"\n")


        for v in range(len(coords_evo_vertex)):
            #print(v)
            NeighR, Neigh_C = fc.find_vertex_neighbour_centers(vorRegions,vorPointRegion,v)
            Neigh_V = fc.find_vertex_neighbour_vertices(vorRidges,v)
            with open('movieoutput/vertices'+str(i)+'.txt','a') as text:
                text.write(str(v)+";"+str(coords_evo_vertex[v,0])+";"+str(coords_evo_vertex[v,1])+";"+str(Neigh_C)+";"+str(Neigh_V)+"\n")
                
        for e in range(len(vorRidges)):
            with open('movieoutput/edges'+str(i)+'.txt','a') as text:
                text.write(str(e)+";"+str(vorRidges[e])+"\n")

        #First find the boundary
        boundary_tissue = fc.find_boundary_vertices(len(coords_evo_vertex),vorRidges)
        truebound = np.where(coords_evo_vertex[boundary_tissue,0]**2+coords_evo_vertex[boundary_tissue,1]**2>9)
        boundary_tissue= list(np.array(boundary_tissue)[truebound])
        boundary_wound = vorRegions[wloc]

        with open('movieoutput/boundaries'+str(i)+'.txt','a') as text:
            text.write(str(boundary_tissue)+"\n")
            text.write(str(boundary_wound)+"\n")


        with open('movieoutput/woundloc'+str(i)+'.txt','a') as text:
            text.write(str(wloc)+"\n")
        
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

    
