import numpy as np
import scipy as scp
import selfpropelledparticlevoronoi as spv
import findconnections as fc
import topologicaltransitions as tpt

def adjacency_matrix(n_vertices, ridges):
    '''
    Calculates the adjacency matrix of the network
    '''
    
    ridges = list(ridges)
    A_ij = np.zeros((n_vertices,n_vertices))
    for ridge in ridges:
        if (ridge[0]!= -1) and (ridge[1]!= -1):
            A_ij[int(ridge[0]),int(ridge[1])] = 1
            A_ij[int(ridge[1]),int(ridge[0])] = 1
        
    return A_ij

def weight_matrix(vertices, ridges):
    '''
    Calculates the weight matrices of the network
    '''
    
    n_vertices = len(vertices)
    
    ridges = list(ridges)
    W_ij = np.zeros((n_vertices,n_vertices))
    for ridge in ridges:
        if (ridge[0]!= -1) and (ridge[1]!= -1):
            W_ij[int(ridge[0]),int(ridge[1])] = fc.norm(vertices[ridge[1]]-vertices[ridge[0]])
            W_ij[int(ridge[1]),int(ridge[0])] = fc.norm(vertices[ridge[0]]-vertices[ridge[1]])
        
    return W_ij

def graph_energy(Amatrix):
    eig_spectrum = np.linalg.eig(Amatrix)[0]
    return np.sum(np.abs(eig_spectrum))