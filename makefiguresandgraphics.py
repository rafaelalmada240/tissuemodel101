import numpy as np
import matplotlib.pyplot as plt
import selfpropelledparticlevoronoi as sppv
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
Force_center_vector = np.array(data_set['forces_center'])
Force_vertex_vector = np.array(data_set['forces_vertices'])

# Do some initial preprocessing

N = coords_evo.shape[0]

vor1 = Voronoi(coords_evo[:,:,0])
vor2 = Voronoi(coords_evo[:,:,-1])

P1 = np.zeros(N)
P2 = np.zeros(N)
for i in range(N):
    P1[i] = sppv.perimeter_vor(vor1,i)
    P2[i] = sppv.perimeter_vor(vor2,i)
    
L1 = np.zeros(N)
L2 = np.zeros(N)
for i in range(N):
    L1[i] = sppv.nsides_vor(vor1,i)
    L2[i] = sppv.nsides_vor(vor2,i)
    
F_center_array = abs_value(Force_center_vector)
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
plt.hist(L1,10,range=(0,10),alpha=0.5)
plt.hist(L2,10,range=(0,10),alpha=0.5)
plt.xlabel('Number of sides')

plt.subplot(133)
plt.hist(np.abs(coords_evo[:,0,0]-coords_evo[:,0,-1]),alpha=0.5)
plt.hist(np.abs(coords_evo[:,1,0]-coords_evo[:,1,-1]),alpha=0.5)
plt.xlabel('Total displacement')
plt.savefig('GeometricDistributions'+namefig+'.png',dpi=150)
plt.close()


# FIGURE 2
namefig2 = input('Name of second figure/ Number of figure: ')
plt.figure(figsize=(16,6))

plt.subplot(121)

plt.plot(np.median(np.abs(np.diff(coords_evo_vertex[:,0],1)),axis=0),'r-')
plt.plot(np.median(np.abs(np.diff(coords_evo[:,0],1)),axis=0),'k-')

plt.plot(np.median(np.abs(np.diff(coords_evo_vertex[:,1],1)),axis=0),'r--')
plt.plot(np.median(np.abs(np.diff(coords_evo[:,1],1)),axis=0),'k--')

plt.xscale('log')
plt.yscale('log')

plt.xlim(1,coords_evo.shape[2]-1)

plt.title('Average absolute displacement $|\\Delta\\vec{r}|$')
plt.legend(['Vertex X','Center X','Vertex Y','Center Y'])
plt.xlabel('Time (steps)')

plt.subplot(122)
plt.plot(np.median(F_center_array,0),'k')

plt.plot(np.median(F_vertex_array,0),'r')

plt.plot(F_center_array[:,1:].T,'k-',alpha=0.01)
plt.plot(F_vertex_array[:,1:].T,'r-',alpha=0.01)
plt.xscale('log')
plt.yscale('log')

#plt.ylim(0.01,1)
plt.xlim(1,coords_evo.shape[2]-1)
plt.title('Average absolute force $|\\vec{F}|$')
plt.legend(['Vertex','Center'])
plt.xlabel('Time (steps)')
plt.savefig('AverageDisplacementsForces'+namefig2+'.png',dpi=150)
plt.close()



plt.show()