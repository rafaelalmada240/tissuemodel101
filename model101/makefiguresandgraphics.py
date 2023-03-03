import numpy as np
import matplotlib.pyplot as plt
import selfpropelledparticlevoronoi as sppv
import topologicaltransitions as tpt
from scipy.spatial import Voronoi, voronoi_plot_2d
import pathlib
import importlib
import h5py as h5

'''Make some useful figures from the dataset resulting from simulations'''

def abs_value(init_array):
    # N*D*M = for a lattice array, N is the number of points, D is the number of spatial dimensions, M is the number of timesteps
    
    return np.sqrt(np.sum(init_array**2,axis=1))


# Load dataset to use

main_path = pathlib.Path().absolute()
datafileloadname = input('Name of savefile: ')
datafileloadname = datafileloadname + '.h5'
data_set = h5.File(str(main_path) + '/' + datafileloadname,'r')

coords_evo = np.array(data_set['centers'])
coords_evo_vertex = np.array(data_set['vertices'])

Force_vertex_vector = np.array(data_set['forces_vertices'])
PerimeterArray = np.array(data_set['perimeter'])
AreaArray = np.array(data_set['area'])

AdjacencyMatrix = np.array(data_set['AdjacencyMatrix'])
# Do some initial preprocessing

P1 = PerimeterArray[:,0]
P2 = PerimeterArray[:,-1]

    
L1 = AreaArray[:,0]
L2 = AreaArray[:,-1]

    

F_vertex_array = abs_value(Force_vertex_vector)
    
# FIGURE 1
namefig = input('Name of first figure/ Number of figure: ')
plt.figure(figsize=(16,6))
plt.suptitle('Relative frequency')

plt.subplot(131)
plt.hist(P1,25,range=(0,10),alpha=0.5)
plt.hist(P2,25,range=(0,10),alpha=0.5)
plt.xlabel('Perimeter')


plt.subplot(132)
plt.hist(L1,10,range=(0,50),alpha=0.5)
plt.hist(L2,10,range=(0,50),alpha=0.5)
plt.xlabel('Area')

plt.subplot(133)
plt.hist(np.abs(coords_evo[:,0,0]-coords_evo[:,0,-1]),alpha=0.5)
plt.hist(np.abs(coords_evo[:,1,0]-coords_evo[:,1,-1]),alpha=0.5)
plt.xlabel('Total displacement')
plt.savefig('GeometricDistributions'+namefig+'.png',dpi=150)
plt.close()


# FIGURE 2
namefig2 = input('Name of second figure/ Number of figure: ')


disp = np.zeros((coords_evo_vertex.shape[0],coords_evo_vertex.shape[2]))

diff_evo = np.diff(coords_evo_vertex[:,:,:],2)
for i in range(coords_evo_vertex.shape[2]-2):
    for j in range(coords_evo_vertex.shape[0]):
        disp[j,i] = sppv.norm(diff_evo[j,:,i])
plt.figure(figsize=(16,6))

plt.subplot(121)

plt.plot(np.median(disp,axis=0),'r-')
plt.plot(disp.T,'r',alpha=0.01)
#plt.plot(np.median(np.abs(np.diff(coords_evo[:,0],1)),axis=0),'k-')

#plt.plot(np.median(np.abs(np.diff(coords_evo_vertex[:,1],1)),axis=0),'r--')
#plt.plot(np.median(np.abs(np.diff(coords_evo[:,1],1)),axis=0),'k--')

#plt.xscale('log')
#plt.yscale('log')

ymin = np.min(np.median(disp,axis=0)[1:coords_evo.shape[2]-2])
ymax = np.max(np.median(disp,axis=0)[1:coords_evo.shape[2]-2])

plt.ylim(ymin*0.99,ymax*1.01)
plt.xlim(1,coords_evo.shape[2]-2)

plt.title('Average absolute displacement $|\\Delta\\vec{r}|$')
plt.legend(['Vertex $\\hat{X}$'])
plt.xlabel('Time (steps)')#Force_center_vector = np.array(data_set['forces_center'])

plt.subplot(122)
plt.plot(np.median(F_vertex_array,0),'r')
plt.plot(F_vertex_array[:,1:-1].T,'r-',alpha=0.01)
#plt.xscale('log')
#plt.yscale('log')

ymin = np.min(np.median(F_vertex_array,0)[1:coords_evo.shape[2]-2])
ymax = np.max(np.median(F_vertex_array,0)[1:coords_evo.shape[2]-2])

plt.ylim(ymin*0.99,ymax*1.01)
plt.xlim(1,coords_evo.shape[2]-2)
plt.title('Average absolute force $|\\vec{F}|$')
plt.legend(['Vertex $\\hat{X}$'])
plt.xlabel('Time (steps)')
plt.savefig('AverageDisplacementsForces'+namefig2+'.png',dpi=150)
plt.close()



plt.show()

# FIGURE 3
eig_spectrum = []
energy = []
for i in range(AdjacencyMatrix.shape[2]):
    eig_spectrum.append(np.linalg.eig(AdjacencyMatrix[:,:,i])[0])
    energy.append(tpt.graph_energy(AdjacencyMatrix[:,:,i]))
 
Nbar = AdjacencyMatrix.shape[0]   
eig_spectrum = np.array(eig_spectrum)
energy = np.array(energy)

namefig = input('Name of first figure/ Number of figure: ')
plt.figure(figsize=(16,6))
plt.suptitle('Spectral properties')

plt.subplot(131)
plt.bar(np.arange(1,Nbar+1),np.mean(np.abs(eig_spectrum[:-1]),0).T,yerr = np.std(np.abs(eig_spectrum[:-1]),0).T)
plt.title('Eigenspectrum $|\lambda|$')

plt.subplot(132)
plt.plot(np.max(np.abs(eig_spectrum[:-1]),axis=1))
plt.xlabel('Time')
plt.title('Spectral radius')

plt.subplot(133)
plt.plot(energy[:-1])
plt.xlabel('Time')
plt.title('Graph energy')
plt.savefig('GraphSpectralProperties'+namefig+'.png',dpi=150)
plt.close()

data_set.close()
