import numpy as np
from vertexmodelpack import sppvertex as sppv
from vertexmodelpack import randomlattice as rlt
from vertexmodelpack import topologicaltransitions as tpt
from vertexmodelpack import readTissueFiles as rTF
import time

# def compute_forces(vorRegions, vorPointRegion, vorRidges, K_run, A0_run, G_run, Lr, Lw, coords_evo_vertex, coords_evo, wloc, h, boundaries):
#     # Implement the force computation using Numba-optimized code
#     F_vertex = sppv.force_vtx_elastic_wound(vorRegions, vorPointRegion, vorRidges, K_run, A0_run,G_run, Lr, Lw, coords_evo_vertex, coords_evo, wloc, h, boundaries)
#     return F_vertex

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
        
    def apply_shear(self, shear_factor):
        self.vertices[:, 0] += shear_factor * self.vertices[:, 1]

class SimulationParameters:
    def __init__(self, epsilonx, epsilont, K_run, G_run, T, p0_min, p0_max, p0_Nbin, lw_min, lw_max, lw_Nbin, UseT1, simple_output):
        self.epsilonx = epsilonx
        self.epsilont = epsilont
        self.K_run = K_run
        self.G_run = G_run
        self.T = T
        self.p0_min = p0_min
        self.p0_max = p0_max
        self.p0_Nbin = p0_Nbin
        self.lw_min = lw_min
        self.lw_max = lw_max
        self.lw_Nbin = lw_Nbin
        self.UseT1 = UseT1
        self.simple_output = simple_output
        self.L_List = list(np.linspace(p0_min, p0_max, p0_Nbin))
        self.Lw_List = list(np.linspace(lw_min, lw_max, lw_Nbin))

class VertexModelSimulation:
    def __init__(self, tissue, sim_params):
        self.tissue = tissue
        self.sim_params = sim_params
        self.setup_simulation()
        
    def setup_simulation(self):
        self.h = self.sim_params.epsilonx * (5 - (-5)) / (2 * np.sqrt(self.tissue.N))
        self.A0_run = self.tissue.av
        self.mu = 1
        self.dt = (self.sim_params.K_run * np.median(self.tissue.av)) / self.mu * self.sim_params.epsilont
        self.M = int(self.sim_params.T / self.dt)
        self.coords_evo = self.tissue.coords
        self.coords_evo_vertex = self.tissue.vertices
        self.areaWound = self.tissue.areaWound0
    
    def run(self):
        for lr in self.sim_params.L_List:
            for lw in self.sim_params.Lw_List:
                self.run_simulation_for_parameters(lr, lw)
                
    def run_simulation_for_parameters(self, lr, lw):
        Lr = lr * self.sim_params.G_run * 2 * self.sim_params.K_run * np.median(self.tissue.av) ** 0.5
        Lw = (lw - lr) * self.sim_params.K_run * np.median(self.tissue.av) ** 1.5
        i = 0
        periWList = []
        areaWList = []
        transitionsList = []
        total_transitions = 0
        
        
        while i < self.M and self.tissue.areaWound0 / 8 <= self.areaWound <= 8 * self.tissue.areaWound0:
            perimeterWound = sppv.perimeter_vor(self.tissue.vorPointRegion, self.tissue.vorRegions, self.coords_evo_vertex, self.tissue.vorRidges, self.tissue.wloc)
            self.areaWound = sppv.area_vor(self.tissue.vorPointRegion, self.tissue.vorRegions, self.coords_evo_vertex, self.tissue.vorRidges, self.tissue.wloc)
            
            J = 0.01
            Rand_vert = J*np.random.rand(self.coords_evo_vertex.shape[0],2)
            
            F_vertex = sppv.force_vtx_elastic_wound(self.tissue.vorRegions, self.tissue.vorPointRegion, self.tissue.vorRidges, self.sim_params.K_run, self.A0_run, self.sim_params.G_run, Lr, Lw, self.coords_evo_vertex, self.coords_evo, self.tissue.wloc, self.h, self.tissue.Boundaries[0])
            A_vertex = rlt.newWhere(self.coords_evo_vertex + self.mu * F_vertex * self.dt, 6)
            self.coords_evo_vertex += self.mu * A_vertex * (F_vertex + Rand_vert)* self.dt
            self.coords_evo = sppv.cells_avg_vtx(self.tissue.vorRegions, self.tissue.vorPointRegion, np.array(self.coords_evo), np.array(self.coords_evo_vertex))
            
            if self.sim_params.UseT1:
                self.tissue.vorRidges, self.coords_evo_vertex, self.tissue.vorRegions, transition_counter = tpt.T1transition2(np.array(self.coords_evo_vertex), np.array(self.tissue.vorRidges), self.tissue.vorRegions, self.tissue.vorPointRegion, 0.01 * self.tissue.r0)
                total_transitions += transition_counter
            
            periWList.append(perimeterWound)
            areaWList.append(self.areaWound)
            transitionsList.append(transition_counter)
            i += 1
        
        self.i_final = i
        self.save_output(lr, lw, periWList, areaWList, transitionsList)

    def save_output(self, lr, lw, periWList, areaWList, transitionsList):
        foldername = self.tissue.foldername
        if self.sim_params.simple_output:
            rTF.simpleOutputTissues(foldername, [self.tissue.areaWound0, self.sim_params.G_run, lr, lw, self.tissue.N], [periWList, areaWList, transitionsList])
        else:
            rTF.movieOutputTissues(foldername, [self.tissue.areaWound0, self.sim_params.G_run, lr, lw, self.tissue.N], [self.coords_evo, self.tissue.vorPointRegion, self.tissue.vorRegions, self.coords_evo_vertex, self.tissue.vorRidges, self.tissue.Boundaries, self.tissue.wloc])

def main():
    list_inputs = []
    with open('inputrwh.txt', 'r') as text:
        for line in text:
            list_inputs.append(line)

    input_tissues = list_inputs[0]
    tissues_list = input_tissues.split()
    tissues = [int(num) for num in tissues_list]
    
    if bool(int(list_inputs[5])):
        sim_params = SimulationParameters(
            epsilonx=float(list_inputs[1]),
            epsilont=float(list_inputs[2]),
            K_run=float(list_inputs[3]),
            G_run=float(list_inputs[4]),
            T=float(list_inputs[8]),
            p0_min=float(list_inputs[13]),
            p0_max=float(list_inputs[14]),
            p0_Nbin=int(list_inputs[15]),
            lw_min=float(list_inputs[16]),
            lw_max=float(list_inputs[17]),
            lw_Nbin=int(list_inputs[18]),
            UseT1=bool(int(list_inputs[9])),
            simple_output=bool(int(list_inputs[10]))
        )

    bigfoldername = list_inputs[12].strip()
    woundsize = int(list_inputs[11])
    
    for tissue_n in tissues:
        tme = time.time()
        print("tissue number: " + str(tissue_n))
        
        foldername = bigfoldername + '/tissue' + str(tissue_n) + '/size' + str(woundsize)
        tissue = Tissue(foldername, tissue_n, woundsize)
        tissue.apply_shear(shear_factor=0.0)  # Apply shear transformation
        
        simulation = VertexModelSimulation(tissue, sim_params)
        simulation.run()
        
        tf = time.time() - tme
        with open(bigfoldername + '/log' + str(tissue_n) + '.txt', 'a') as log:
            log.write('Median cell area - ')
            log.write(' ')
            log.write(str(np.median(tissue.av).round(6)) + '\n')
            log.write('Spatial resolution of the tissue - ')
            log.write(' ')
            log.write(str(simulation.h.round(6)) + '\n')
            log.write('Time elapsed: ')
            log.write(str(tf) + '\n')
            log.write('Total number of iterations - ')
            log.write(' ')
            log.write(str(simulation.i_final)+'\n')

if __name__ == '__main__':
    main()
