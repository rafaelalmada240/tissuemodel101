import numpy as np
from vertexmodelpack import sppvertex as sppv
from vertexmodelpack import randomlattice as rlt
from vertexmodelpack import topologicaltransitions as tpt
from vertexmodelpack import readTissueFiles as rTF
import time

class Tissue:
    def __init__(self, foldername, tissue_n, woundsize):
        self.foldername = foldername
        self.tissue_n = tissue_n
        self.woundsize = woundsize
        self.load_tissue()
    
    def load_tissue(self):
        dataset = rTF.open_tissuefile(self.foldername, 0)
        self.coords = dataset['centers']
        self.vorPointRegion = dataset['point regions']
        self.vorRegions = dataset['regions']
        self.vertices = dataset['vertices']
        self.vorRidges = dataset['Edge connections']
        self.Boundaries = dataset['boundaries']
        self.wloc = dataset['WoundLoc']
        self.N = len(self.coords[:, 0])
        self.av = sppv.areas_vor(self.vorPointRegion, self.vorRegions, self.vertices, self.vorRidges, self.vorPointRegion)
        self.r0 = np.sqrt(np.median(self.av) / np.pi)
        self.areaWound0 = sppv.area_vor(self.vorPointRegion, self.vorRegions, self.vertices, self.vorRidges, self.wloc)
        
        
class Substrate:
    def __init__(self):
        self.Substrate = []


list_inputs = []
with open('inputrwh.txt','r') as text:
    for line in text:
        list_inputs.append(line)

input_tissues = list_inputs[0]
tissues_list = input_tissues.split()
tissues = [int(num) for num in tissues_list]#[1,2,3,4,5,6,7,8,9,10]


epsilonx = float(list_inputs[1])
epsilont = float(list_inputs[2])

K_run = float(list_inputs[3])
G_run = float(list_inputs[4])

list_plw = bool(int(list_inputs[5]))
if list_plw:
    p0_min = float(list_inputs[13])
    p0_max = float(list_inputs[14])
    p0_Nbin = int(list_inputs[15])

    lw_min = float(list_inputs[16])
    lw_max = float(list_inputs[17])
    lw_Nbin = int(list_inputs[18])

    #To make a list of Ls and Lws (dimensionless)
    L_List =  list(np.linspace(p0_min,p0_max,p0_Nbin))
    print("bin resolution p0")
    print(float((p0_max-p0_min)/(p0_Nbin-1)))


    Lw_List = list(np.linspace(lw_min,lw_max,lw_Nbin))
    print("bin resolution lw")
    print(float((lw_max-lw_min)/(lw_Nbin-1)))
    
else:
    p0_min = float(list_inputs[6])
    p0_max = p0_min
    p0_Nbin = 2

    lw_min = float(list_inputs[7])
    lw_max = lw_min
    lw_Nbin = 2

    #To make a list of Ls and Lws (dimensionless)
    L_List =  [p0_min]
    print("bin resolution p0")
    print(float((p0_max-p0_min)/(p0_Nbin-1)))


    Lw_List = [lw_min]
    print("bin resolution lw")
    print(float((lw_max-lw_min)/(lw_Nbin-1)))
    

L_max = 5
L_min = -5  
DeltaL = L_max - L_min

T = float(list_inputs[8])

UseT1 = bool(int(list_inputs[9]))

simple_output = bool(int(list_inputs[10]))

woundsize = int(list_inputs[11]) 

# Opening files in respective folders - probably use a function that turns this into a dictionary:
##########################################################################################################
bigfoldername = list_inputs[12].strip()

for tissue_n in tissues:
    
    tme = time.time()
    
    print("tissue number: "+str(tissue_n))
    
    foldername = bigfoldername + '/tissue' + str(tissue_n) + '/size' + str(woundsize)
    tissue = Tissue(foldername, tissue_n, woundsize)
    print(tissue.N)

#################################################################################################################
    
    
    # Simulation parameters    
    # av = []
    
    # av = sppv.areas_vor(vorPointRegion,vorRegions,vertices,vorRidges,vorPointRegion)
    print("Median cell area")
    print(np.median(tissue.av))
    # print(np.std(av))
   
    # r0 = np.sqrt(np.median(av)/np.pi)
    h = epsilonx*DeltaL/(2*np.sqrt(tissue.N)) #size step
    print("Simulation spatial resolution")
    print(h)
    

    A0_run = tissue.av# or np.median(av), maybe it depends on the model
    mu = 10
    dt = (K_run*np.median(tissue.av))/mu*epsilont #time step normalized by K_run A0_run/mu
    print("Simulation temporal resolution")
    print(dt)

    #for PC
    
    #for cluster
    #M = 50000 
    

    M = int(T/dt)
    print("Simulation Max number of iterations")
    print(M)
    


#####################################################################################################################
    i_final = 0
    for lr in L_List:
        # Run simulations
        for lw in Lw_List:
            Lr = lr*G_run*2*K_run*np.median(tissue.av)**(1/2)
            Lw = (lw-lr)*K_run*np.median(tissue.av)**(3/2)
            
            i = 0
            transitionsList = []
            periWList = []
            areaWList= []
            total_transitions = 0
            areaWound = tissue.areaWound0
            coords_evo = tissue.coords
            coords_evo_vertex = tissue.vertices
            
            list_coords = []
            list_vertex = []
            list_points = []
            list_regions = []
            list_edges = []
            list_boundaries = []
            list_wloc = []

            while (i < M) and ((areaWound>= tissue.areaWound0/16) and (areaWound<= 16*tissue.areaWound0)):
                
                #Compute the area and perimeter of the wound
                perimeterWound = sppv.perimeter_vor(tissue.vorPointRegion,tissue.vorRegions,coords_evo_vertex,tissue.vorRidges,tissue.wloc)
                areaWound = sppv.area_vor(tissue.vorPointRegion,tissue.vorRegions,coords_evo_vertex,tissue.vorRidges,tissue.wloc)
                
                #Compute forces
                J = 0.5
                
                Rand_vertex = J*(np.random.rand(coords_evo_vertex.shape[0],2)-0.5)
                F_vertex = sppv.force_vtx_elastic_wound(tissue.vorRegions, tissue.vorPointRegion, tissue.vorRidges, K_run,A0_run,G_run,Lr,Lw,coords_evo_vertex,coords_evo,tissue.wloc,h,tissue.Boundaries[0])
                # print(np.mean(F_vertex,0))    
                #Reflexive boundary conditions
                A_vertex = rlt.newWhere(coords_evo_vertex + mu*F_vertex*dt,30)
                coords_evo_vertex = coords_evo_vertex + mu*A_vertex*(F_vertex+Rand_vertex)*dt
                
                # if i == 10:
                #     print(F_vertex[int(np.random.rand()*200)])
                #     print(Rand_vertex[int(np.random.rand()*200)])
                #     print(F_vertex[tissue.vorRegions[tissue.wloc][0]])
                #     print(Rand_vertex[tissue.vorRegions[tissue.wloc][0]])
                #     tissue.wloc[0]
                
            
                #Cell center positions are the average of the cell vertex positions
                coords_evo = sppv.cells_avg_vtx(tissue.vorRegions,tissue.vorPointRegion,np.array(coords_evo),np.array(coords_evo_vertex))
                
                #Do topological rearrangements
                transition_counter = 0
                if UseT1:
                    tissue.vorRidges, coords_evo_vertex, tissue.vorRegions, transition_counter =  tpt.T1transition2(np.array(coords_evo_vertex),np.array(tissue.vorRidges),tissue.vorRegions, tissue.vorPointRegion,0.1*tissue.r0)
                    
                i = i + 1
                total_transitions += transition_counter 
                
                #Store values in list to be saved later
                periWList.append(perimeterWound)
                areaWList.append(areaWound)
                transitionsList.append(transition_counter)
                
                list_coords.append(coords_evo)
                list_vertex.append(coords_evo_vertex)
                list_edges.append(tissue.vorRidges)
                list_points.append(tissue.vorPointRegion)
                list_regions.append(tissue.vorRegions)
                list_boundaries.append(tissue.Boundaries)
                list_wloc.append(tissue.wloc)
            
            
            #Output of simulations for analysis
            
            if simple_output: 
                rTF.simpleOutputTissues(foldername,[tissue.areaWound0,G_run,lr,lw,tissue.N],[periWList,areaWList,transitionsList])
            else:
                rTF.movieOutputTissues(foldername,[len(list_wloc),lr,lw],[list_coords,list_points,list_regions,list_vertex,list_edges,list_boundaries,list_wloc])
                
            i_final += i
    #################################################################################################################################################################   
            
    #Log of simulations

    tf = time.time()-tme

    with open(bigfoldername+'/log'+str(tissue_n)+'.txt','a') as log:
        log.write('Median cell area - ')
        log.write(' ')
        log.write(str(np.median(tissue.av).round(6))+'\n')
        log.write('Spatial resolution of the tissue - ')
        log.write(' ')
        log.write(str(h.round(6))+'\n')
        log.write('Temporal resolution of the tissue - ')
        log.write(' ')
        log.write(str(dt.round(6))+'\n')
        log.write('Simulation time (s) - ')
        log.write(' ')
        log.write(str(tf)+'\n')
        log.write('Total number of iterations - ')
        log.write(' ')
        log.write(str(i_final)+'\n')
        log.write('Bin resolution p0 - ')
        log.write(' ')
        log.write(str(float((p0_max-p0_min)/(p0_Nbin-1)))+'\n')
        log.write('Bin resolution lw - ')
        log.write(' ')
        log.write(str(float((lw_max-lw_min)/(lw_Nbin-1)))+'\n')
        log.write('Used T1? - ')
        log.write(' ')
        log.write(str(UseT1)+'\n')
