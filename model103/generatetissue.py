import numpy as np
import randomlatticethermalized as rlt
import selfpropelledparticlevoronoi as spv
import h5py as h5
import pathlib

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

if (size_of_wound == 1) or (size_of_wound == 0):
    print('Single cell wound')
else: 
    wound_radius = size_of_wound*dL/2
    L_wound = np.zeros(N)
    for i in range(N):
        L_wound[i] = np.linalg.norm(coords[i]-coords[wound_loc])
    wound_cells_loc = list(np.where(L_wound < wound_radius)[0])
    wound_cells_loc.remove(wound_loc)
    
    coords = np.delete(coords,np.array(wound_cells_loc),0)
    N_new = len(coords)
    L_array = np.zeros(N_new)
    for i in range(N_new):
        L_array[i] = np.linalg.norm(coords[i])
    
    wound_loc = np.argmin(L_array)
    N = N_new
    
    



    

main_path = pathlib.Path().absolute()
datafilesavename = input('Name of savefile: ')
datafilesavename = datafilesavename + '.h5'
data_save = h5.File(datafilesavename,'w')
data_save.create_dataset('centers',data = coords)
data_save.create_dataset('number of centers',data=N) 
data_save.create_dataset('woundlocation',data=wound_loc) 
data_save.close()