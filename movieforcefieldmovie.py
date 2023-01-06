import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
import movies_from_plots as mfp
from scipy import ndimage as nd 
import h5py as h5
import importlib
import pathlib

''' 
Makes a short movie showing the evolution of the network force field over time
'''

def frame_force_field_evo(center_array, force_array):
    
    # Make a figure for each timestep of our simulation
    
    xnew = np.arange(-5.01, 5.01, 5e-2)
    ynew = np.arange(-5.01, 5.01, 5e-2)
    Xnew,Ynew = np.meshgrid(xnew,ynew)
    
    for i in range(len(center_array[0,0])):
        f1 = LinearNDInterpolator(list(zip(center_array[:,0,i],center_array[:,1,i])),force_array[:,i])
        Pnew1 = f1(Xnew, Ynew)
        fig = plt.figure(figsize=(16,6))
    
        plt.subplot(121)
        plt.pcolormesh(Xnew,Ynew,np.log10(nd.gaussian_filter(Pnew1,1)),shading='auto',cmap='seismic')
        plt.plot(center_array[:,0,i],center_array[:,1,i],'g*',markersize = 6)
        
        plt.colorbar()
        plt.clim(-2,0)
        plt.xlabel('x',fontsize = 16)
        plt.ylabel('y',fontsize = 16)
        plt.title('Step = '+str(i),fontsize = 16)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
    
        ax = fig.add_subplot(122,projection='3d')
        ax.plot_surface(Xnew,Ynew,np.log10(nd.gaussian_filter(Pnew1,1)),lw=0.5,rstride = 8, cstride=8,alpha=0.5)
        ax.set(xlim=(-5,5),ylim=(-5,5),zlim=(-2,2),title='Force field (log scale)',xlabel='X',ylabel='Y',zlabel='Z')

        plt.savefig('f'+str(int(i))+'.png',dpi=100)
        plt.close(fig)
    plt.show()    
    return

def abs_value(init_array):
    # N*D*M = for a lattice array, N is the number of points, D is the number of spatial dimensions, M is the number of timesteps
    return np.sqrt(np.sum(init_array**2,axis=1))

# Load a dataset and take the necessary arrays
main_path = pathlib.Path().absolute()
datafileloadname = input('Name of savefile: ')
datafileloadname = datafileloadname + '.h5'
data_set = h5.File(str(main_path) + '/' + datafileloadname,'r')
coords_evo = np.array(data_set['centers'])
Force_center_vector = np.array(data_set['forces_center'])

F_center_array = abs_value(Force_center_vector)

#Make a movie

N_initial = 0 #int(input('Starting frame: '))
N_final = len(coords_evo[0,0]) #int(input('Stopping frame: '))
fr_step = 1 #int(input('Frame step: '))
image_file = 'f'
filename = input('Video file name (with .mp4 included): ')

frame_force_field_evo(coords_evo,F_center_array)
img_array,size = mfp.loadframe(N_initial,N_final,fr_step,image_file)
mfp.savevideo(img_array,filename,size,image_file)