import numpy as np
import findconnections as fc
import networkx as nx
from time import time



'''
    regions - (list) set of all the vertices that compose the different regions of the network
    Some regions are empty, but it is best not to remove them for consistency
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    cells - array of coordinates for all cell centers of the network
    vertices - array of coordinate positions for each vertex of the network

'''

## Functions to calculate perimeter and area terms for a given vertex

def adj_mat(R,ridges):
    arrayR = np.array(R)
    bin_mat = np.zeros((len(R),len(R)))
    for vi in R:
        N_v = fc.find_vertex_neighbour_vertices(ridges,vi)
        loc_i = np.argwhere(arrayR==vi)[0][0]
        for vj in N_v:
            if vj in R:
                loc_j = np.argwhere(arrayR==vj)[0][0]
                bin_mat[loc_i,loc_j] = 1
    
    return bin_mat
    
def flatten(l):
    return [item for sublist in l for item in sublist]

def rearrange(n, bin_mat,to_print = False):
    

    
    G = nx.from_numpy_array(bin_mat)
    cyclesG = nx.cycle_basis(G,0)

    #print("sorted"+str(sorted(cyclesG)))
    if len(cyclesG)>0:
        arr_new = sorted(cyclesG)[0]
        if to_print == True:
            print(arr_new)
    else:
        arr_new = list(np.arange(n))

      
    return arr_new



def perimeter_vor(point_region,regions,vertices,ridges,i):

    R = fc.find_center_region(regions,point_region,i)
    
    #Calculate the perimeter term
    if len(R)>2:
        rearrange_loc = rearrange(len(R),adj_mat(R,ridges))
        if len(rearrange_loc)==len(R):
            V = vertices[R][rearrange_loc]
        else:
            V = vertices[R]
        P = 0
        for i in range(len(V)):
            P += fc.norm(V[(i+1)%len(V)]-V[i])
    else:
        P = 0
    return P


def area_vor(point_region,regions,vertices,ridges,i,to_print = False):

    R = fc.find_center_region(regions,point_region,i)

    
    #Calculate the area term
    if len(R)>2:
        rearrange_loc = rearrange(len(R),adj_mat(R,ridges),to_print)
    
        if len(rearrange_loc)==len(R):
            V = vertices[R][rearrange_loc]
        else:
            V = vertices[R]
        A1 = 0
        for i in range(len(V)):
            A1 += np.cross(V[i],V[(i+1)%len(V)])  
    else:
        A1 = 0
                
    return 1/2*fc.norm(A1)

def perimeters_vor(point_region,regions,vertices,ridges,list_i):

    Rlist = [fc.find_center_region(regions,point_region,i) for i in list_i]
    Plist = []
    for R in Rlist:
    #Calculate the perimeter term
        if len(R)>2:
            rearrange_loc = rearrange(len(R),adj_mat(R,ridges))
            if len(rearrange_loc)==len(R):
                V = vertices[R][rearrange_loc]
            else:
                V = vertices[R]
            P = np.sum(np.array([fc.norm(V[(i+1)%len(V)]-V[i]) for i in range(len(V))]))
            #for i in range(len(V)):
                #P += 
        else:
            P = 0
        Plist.append(P)
    return Plist

def areas_vor(point_region,regions,vertices,ridges,list_i,to_print = False):
    Rlist = [fc.find_center_region(regions,point_region,i) for i in list_i]
    Alist = []
    for R in Rlist:
        #Calculate the area term
        if len(R)>2:
            rearrange_loc = rearrange(len(R),adj_mat(R,ridges),to_print)
        
            if len(rearrange_loc)==len(R):
                V = vertices[R][rearrange_loc]
            else:
                V = vertices[R]
            A1 = np.sum(np.array([np.cross(V[i],V[(i+1)%len(V)])  for i in range(len(V))]))
            #for i in range(len(V)):
            #    A1 += np.cross(V[i],V[(i+1)%len(V)])  
        else:
            A1 = 0
        Alist.append(1/2*fc.norm(A1))
    return Alist
            
        
def nsides_vor(point_region,regions,i):

    R = fc.find_center_region(regions,point_region,i)
    nsides=len(R)
    return nsides



# def energy_vor_v2(point_region,regions,ridges, vertices,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound):
    
#     R,N_c = fc.find_vertex_neighbour_centers(regions, point_region,vertex)
#     N_v = fc.find_vertex_neighbour_vertices(ridges,vertex)
    
#     Intersect_nv = list(set(N_v).intersection(bound_wound))
#     E = 0
#     # tm = time() 
#     # for i in N_c:
#     #     if i != wloc:
#     #         Pi = perimeter_vor(point_region,regions,vertices,ridges, i)
#     #         Ai = area_vor(point_region,regions,vertices, ridges, i)
            
#     #         E += K/2*(Ai-A0[i])**2 + G/2*Pi**2
        
        
#     Ncw = list(N_c)
#     if wloc in Ncw:
#         Ncw.remove(wloc)
       
        
#     P = perimeters_vor(point_region,regions,vertices,ridges, Ncw)#[perimeter_vor(point_region,regions,vertices,ridges, i) for i in N_c if i != wloc]
#     A = areas_vor(point_region,regions,vertices, ridges, Ncw)#[area_vor(point_region,regions,vertices, ridges, i) for i in N_c if i != wloc]
#     # print(time()-tm)
#     A0c = [A0[i] for i in Ncw]
#     ESum = np.array([K/2*(A[i]-A0c[i])**2 + G/2*P[i]**2 for i in range(len(A))])
    
    
#     E = E + np.sum(ESum)
#     for j in N_v:
#         if (vertex not in boundary_tissue) and (j not in boundary_tissue):
#             v = vertices[j]        
#             edgeV = vertices[vertex] - v
#             lj = fc.norm(edgeV)    
#             if j not in Intersect_nv:    
#                 E += L*lj
#             else:
#                 E += -Lw*lj
    
    
#     return E

# def force_vtx_finite_gradv2(point_region, regions, ridges, vertices, vertex, K, A0, G, L,h,boundary_tissue, Lw,wloc, bound_wound):
    
#     new_vertices1x = np.array(vertices)
#     new_vertices2x = np.array(vertices)
        
#     new_vertices1y = np.array(vertices)
#     new_vertices2y = np.array(vertices)
        
#     f_v = 0.0*np.array([1.,1.])
        
#     n1 = np.array([1,0])
#     n2 = np.array([0,1])
        
#     new_vertices1x[vertex] = new_vertices1x[vertex] - h*n1
#     new_vertices2x[vertex] = new_vertices2x[vertex] + h*n1
        
#     new_vertices1y[vertex] = new_vertices1y[vertex] - h*n2
#     new_vertices2y[vertex] = new_vertices2y[vertex] + h*n2
        
#     Ev1x = energy_vor_v2(point_region,regions,ridges, new_vertices1x,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
#     Ev2x = energy_vor_v2(point_region,regions,ridges, new_vertices2x,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
        
#     Ev1y = energy_vor_v2(point_region,regions,ridges, new_vertices1y,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
#     Ev2y = energy_vor_v2(point_region,regions,ridges, new_vertices2y,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
        
#     dEdx = 0.5*(Ev2x-Ev1x)/h
#     dEdy = 0.5*(Ev2y-Ev1y)/h
        
#     f_v = -(dEdx*n1 + dEdy*n2)
    
#     return f_v


# def force_vtx_elastic_wound(regions,point_region, ridges, K,A0,G,L,Lw,vertices,centers, wloc,h):
    
#     '''
    
#     Calculates the force in all of the vertices according to the vertex model (energy gradient descent method)
#     Accounting for a wound healing model
    
#     Variables:
    
#     regions - (list) set of all the vertices that compose the different regions of the network
#     point_region - (list) set of the different regions of the network
#     ridges - (list) set of all edges in the network
#     K, A0, G, L - (float) model parameters 
#     vertices - array of coordinate positions for each vertex of the network
#     h - step of gradient descent
    
#     Output:

#     F_V - array of all forces acting on the vertices

    
#     '''
    
#     LV = len(vertices) #Number of vertices

#     F_V = [] # Vector of resulting forces acting on the vertices
    
    
#     r0 = np.sqrt(np.mean(A0)/np.pi)
    
#     #For vertices, well, this may be a huge mess
    
#     #First find the boundary
#     boundary_tissue = fc.find_boundary_vertices(len(vertices),ridges)
#     truebound = np.where(vertices[boundary_tissue,0]**2+vertices[boundary_tissue,1]**2>9)
#     boundary_tissue= np.array(boundary_tissue)[truebound]
#     boundary_wound = fc.find_wound_boundary(regions,point_region,wloc)

#     for v in range(LV):
#         if v not in boundary_tissue:
#             f_v = force_vtx_finite_gradv2(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue,Lw,wloc, boundary_wound)
#             #else:
#             #    f_v = force_vtx_finite_grad_wound(point_region, regions, ridges, vertices, v, K, A0, G, L,)
#                 #print(fc.norm(f_v))
#             NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
#             for i in range(len(NeighC)):
#                 for j in range(i+1,len(NeighC)):
#                     ci = centers[NeighC[i]]
#                     cj = centers[NeighC[j]]
#                     rij = fc.norm(cj-ci)
#                     nij = (cj-ci)/(rij+h)
#                     if rij <= r0:
#                         f_v += 0.1*(rij-r0)*nij
#                     else:
#                         continue
#         else:
#             f_v = 0.0*np.array([1.,1.]) 
                
#         #Maybe include a regularizing force acting on the cells
#         F_V.append(f_v)
        
#     return np.array(F_V)


# def cells_avg_vtx(regions,point_region,cells,vertices):
#     for i in range(len(cells)):
#         Neigh_c = fc.find_center_region(regions,point_region, i)
#         avg_vc = np.mean(vertices[Neigh_c],0)
#         cells[i] = avg_vc
        
#     return cells


#########################################################################################################################################3

def energy_vor_v2(point_region,regions,ridges, vertices,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound):
    
    R,N_c = fc.find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = fc.find_vertex_neighbour_vertices(ridges,vertex)
    
    Intersect_nv = list(set(N_v).intersection(bound_wound))
    E = 0
    for i in N_c:
        if i != wloc:
            Pi = perimeter_vor(point_region,regions,vertices,ridges, i)
            Ai = area_vor(point_region,regions,vertices, ridges, i)
            
            E += K/2*(Ai-A0[i])**2 + G/2*Pi**2
    for j in N_v:
        if (vertex not in boundary_tissue) and (j not in boundary_tissue):
            v = vertices[j]        
            edgeV = vertices[vertex] - v
            lj = fc.norm(edgeV)    
            if j not in Intersect_nv:    
                E += L*lj
            else:
                E += -Lw*lj
    
    
    return E

def force_vtx_finite_gradv2(point_region, regions, ridges, vertices, vertex, K, A0, G, L,h,boundary_tissue, Lw,wloc, bound_wound):
    
    new_vertices1x = np.array(vertices)
    new_vertices2x = np.array(vertices)
        
    new_vertices1y = np.array(vertices)
    new_vertices2y = np.array(vertices)
        
    f_v = 0.0*np.array([1.,1.])
        
    n1 = np.array([1,0])
    n2 = np.array([0,1])
        
    new_vertices1x[vertex] = new_vertices1x[vertex] - h*n1
    new_vertices2x[vertex] = new_vertices2x[vertex] + h*n1
        
    new_vertices1y[vertex] = new_vertices1y[vertex] - h*n2
    new_vertices2y[vertex] = new_vertices2y[vertex] + h*n2
        
    Ev1x = energy_vor_v2(point_region,regions,ridges, new_vertices1x,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
    Ev2x = energy_vor_v2(point_region,regions,ridges, new_vertices2x,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
        
    Ev1y = energy_vor_v2(point_region,regions,ridges, new_vertices1y,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
    Ev2y = energy_vor_v2(point_region,regions,ridges, new_vertices2y,vertex, K,A0,G,L,boundary_tissue,Lw,wloc, bound_wound)
        
    dEdx = 0.5*(Ev2x-Ev1x)/h
    dEdy = 0.5*(Ev2y-Ev1y)/h
        
    f_v = -(dEdx*n1 + dEdy*n2)
    
    return f_v


def force_vtx_elastic_wound(regions,point_region, ridges, K,A0,G,L,Lw,vertices,centers, wloc,h):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (energy gradient descent method)
    Accounting for a wound healing model
    
    Variables:
    
    regions - (list) set of all the vertices that compose the different regions of the network
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    K, A0, G, L - (float) model parameters 
    vertices - array of coordinate positions for each vertex of the network
    h - step of gradient descent
    
    Output:

    F_V - array of all forces acting on the vertices

    
    '''
    
    LV = len(vertices) #Number of vertices

    F_V = [] # Vector of resulting forces acting on the vertices
    
    
    r0 = np.sqrt(np.mean(A0)/np.pi)
    
    #For vertices, well, this may be a huge mess
    
    #First find the boundary
    boundary_tissue = fc.find_boundary_vertices(len(vertices),ridges)
    truebound = np.where(vertices[boundary_tissue,0]**2+vertices[boundary_tissue,1]**2>9)
    boundary_tissue= np.array(boundary_tissue)[truebound]
    boundary_wound = fc.find_wound_boundary(regions,point_region,wloc)

    for v in range(LV):
        if v not in boundary_tissue:
            f_v = force_vtx_finite_gradv2(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue,Lw,wloc, boundary_wound)
            #else:
            #    f_v = force_vtx_finite_grad_wound(point_region, regions, ridges, vertices, v, K, A0, G, L,)
                #print(fc.norm(f_v))
            NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
            for i in range(len(NeighC)):
                for j in range(i+1,len(NeighC)):
                    ci = centers[NeighC[i]]
                    cj = centers[NeighC[j]]
                    rij = fc.norm(cj-ci)
                    nij = (cj-ci)/(rij+h)
                    if rij <= r0:
                        f_v += 0.1*(rij-r0)*nij
                    else:
                        continue
        else:
            f_v = 0.0*np.array([1.,1.]) 
                
        #Maybe include a regularizing force acting on the cells
        F_V.append(f_v)
        
    return np.array(F_V)


def cells_avg_vtx(regions,point_region,cells,vertices):
    for i in range(len(cells)):
        Neigh_c = fc.find_center_region(regions,point_region, i)
        avg_vc = np.mean(vertices[Neigh_c],0)
        cells[i] = avg_vc
        
    return cells