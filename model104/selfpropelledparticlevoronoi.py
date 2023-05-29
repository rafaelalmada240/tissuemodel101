import numpy as np
import findconnections as fc

'''
    regions - (list) set of all the vertices that compose the different regions of the network
    Some regions are empty, but it is best not to remove them for consistency
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    cells - array of coordinates for all cell centers of the network
    vertices - array of coordinate positions for each vertex of the network

'''

## Functions to calculate perimeter and area terms for a given vertex

def perimeter_vor(point_region,regions,vertices,i):

    R = fc.find_center_region(regions,point_region,i)
    V = vertices[R]
    P = 0
    #Calculate the perimeter term
    if len(R)>2:
        for i in range(len(V)):
            P += fc.norm(V[(i+1)%len(V)]-V[i])
    return P

def area_vor(point_region,regions,vertices,i):

    R = fc.find_center_region(regions,point_region,i)
    V = vertices[R]
    A = 0
    #Calculate the area term
    if len(R)>2:
        for i in range(len(V)):
            A += 1/2*fc.norm((np.cross(V[i],V[(i+1)%len(V)])))        
    return A
        
def nsides_vor(point_region,regions,i):

    R = fc.find_center_region(regions,point_region,i)
    nsides=len(R)
    return nsides

def energy_vor(point_region,regions,ridges, vertices,vertex, K,A0,G,L,boundary_tissue):

    R,N_c = fc.find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = fc.find_vertex_neighbour_vertices(ridges,vertex)
    E = 0
    
    #Sum over all cells containing the vertex
    for i in N_c:
        Pi = perimeter_vor(point_region,regions,vertices, i)
        Ai = area_vor(point_region,regions,vertices, i)
        E += K/2*(Ai-A0[i])**2 + G/2*Pi**2
        
    #Sum over all adjacent vertices (if both vertices are not in boundary)
    for j in N_v:
        if (vertex not in boundary_tissue) and (j not in boundary_tissue):
            v = vertices[j]        
            edgeV = vertices[vertex] - v
            lj = fc.norm(edgeV)    
            E += L*lj
        
    return E

def energy_vor_wound(point_region,regions,ridges, vertices,vertex, K,A0,G,L,Lw,wloc, bound_wound):

    R,N_c = fc.find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = fc.find_vertex_neighbour_vertices(ridges,vertex)
    
    Intersect_nv = list(set(N_v).intersection(bound_wound))
    E = 0
    
    #Sum over all cells containing the vertex
    for i in N_c:
        if i != wloc:
            Pi = perimeter_vor(point_region,regions,vertices, i)
            Ai = area_vor(point_region,regions,vertices, i)
            E += K/2*(Ai-A0[i])**2 + G/2*Pi**2
        
    #Sum over all adjacent vertices
    for j in N_v:
        
        v = vertices[j]        
        edgeV = vertices[vertex] - v
        lj = fc.norm(edgeV)
        if j not in Intersect_nv:    
            E += L*lj
        else:
            E += -Lw*lj
        
    return E


def force_vtx(point_region, regions, ridges, vertices, vertex, K, A0, G, L):
    
    R,N_c = fc.find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = fc.find_vertex_neighbour_vertices(ridges,vertex)
    F = 0
    z = np.array([0,0,1])
    
    for i in range(len(N_c)):
        c_i = N_c[i]
        R_i = regions[R[i]]
        if len(R_i)>2:
            S = list(set(R_i).intersection(N_v))
        
            if len(S)>1:
        
                n_a = -np.cross(z,vertices[S[0]]-vertices[S[1]])[:-1]
        
                v_p0 = vertices[vertex]-vertices[S[0]]
                v_p1 = vertices[S[1]]-vertices[vertex]
                l_0 = fc.norm(v_p0)
                l_1 = fc.norm(v_p1)
                n_p = v_p0/l_0 - v_p1/l_1
        
                Pi = perimeter_vor(point_region,regions,vertices, c_i)
                Ai = area_vor(point_region,regions,vertices, c_i)
                
                print(Pi)
                print(Ai)
                F += K/2*(Ai-A0)*n_a - float(G*Pi+L/2)*n_p
        
    return F

def force_vtx_finite_grad(point_region, regions, ridges, vertices, vertex, K, A0, G, L,h,boundary_tissue):
    
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
        
    Ev1x = energy_vor(point_region,regions,ridges, new_vertices1x,vertex, K,A0,G,L,boundary_tissue)
    Ev2x = energy_vor(point_region,regions,ridges, new_vertices2x,vertex, K,A0,G,L,boundary_tissue)
        
    Ev1y = energy_vor(point_region,regions,ridges, new_vertices1y,vertex, K,A0,G,L,boundary_tissue)
    Ev2y = energy_vor(point_region,regions,ridges, new_vertices2y,vertex, K,A0,G,L,boundary_tissue)
        
    dEdx = 0.5*(Ev2x-Ev1x)/h
    dEdy = 0.5*(Ev2y-Ev1y)/h
        
    f_v = -(dEdx*n1 + dEdy*n2)
        
    return f_v


def force_vtx_finite_grad_wound(point_region, regions, ridges, vertices, vertex, K, A0, G, L, Lw,h,wloc, bound_wound):
    
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
        
    Ev1x = energy_vor_wound(point_region,regions,ridges, new_vertices1x,vertex, K,A0,G,L,Lw,wloc, bound_wound)
    Ev2x = energy_vor_wound(point_region,regions,ridges, new_vertices2x,vertex, K,A0,G,L,Lw,wloc, bound_wound)
        
    Ev1y = energy_vor_wound(point_region,regions,ridges, new_vertices1y,vertex, K,A0,G,L,Lw,wloc, bound_wound)
    Ev2y = energy_vor_wound(point_region,regions,ridges, new_vertices2y,vertex, K,A0,G,L,Lw,wloc, bound_wound)
        
    dEdx = 0.5*(Ev2x-Ev1x)/h
    dEdy = 0.5*(Ev2y-Ev1y)/h
        
    f_v = -(dEdx*n1 + dEdy*n2)
        
    return f_v
## Different force models 



     

def force_vtx_elastic(regions,point_region, ridges, K,A0,G,L, vertices,centers, h):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (energy gradient descent method)
    
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
    boundary_tissue = fc.find_boundary_vertices(len(vertices),ridges)

    for v in range(LV):
        if v not in boundary_tissue:
        
            f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
        
            #Maybe include a regularizing force acting on the cells
        
        
            NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
            for i in range(len(NeighC)):
                for j in range(i+1,len(NeighC)):
                    ci = centers[NeighC[i]]
                    cj = centers[NeighC[j]]
                    rij = fc.norm(cj-ci)
                    nij = (cj-ci)/rij
                    if rij <= r0:
                        f_v += (rij-r0)*nij
                    else:
                        continue
        else:
            f_v = 0.0*np.array([1.,1.]) 
        
        F_V.append(f_v)
        
    return np.array(F_V)



def force_vtx_elastic_shear(regions,point_region, ridges, K,A0,G,L, vertices,centers, magF, h, boundary_tissue, bound_up, bound_down):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (energy gradient descent method)
    
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
    #boundary_tissue = fc.find_boundary_vertices(len(vertices),ridges)

    for v in range(LV):
        if v not in boundary_tissue:
        
            f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
        
            #Maybe include a regularizing force acting on the cells
        
        
            NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
            for i in range(len(NeighC)):
                for j in range(i+1,len(NeighC)):
                    ci = centers[NeighC[i]]
                    cj = centers[NeighC[j]]
                    rij = fc.norm(cj-ci)
                    nij = (cj-ci)/rij
                    if rij <= r0:
                        f_v += 0.1*(rij-r0)*nij
                    else:
                        continue
        else:
            if v in bound_down:
                f_v = 0.0*np.array([1.,1.]) 
            else:
                if v in bound_up: 
                    n1 = np.array([1,0])
                    f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)*n1 + magF*n1
                    
                    NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
                    for i in range(len(NeighC)):
                        for j in range(i+1,len(NeighC)):
                            ci = centers[NeighC[i]]
                            cj = centers[NeighC[j]]
                            rij = fc.norm(cj-ci)
                            nij = (cj-ci)/rij
                            if rij <= r0:
                                f_v += 0.1*(rij-r0)*nij
                            else:
                                continue
                else:
                    f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
                    
                    NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
                    for i in range(len(NeighC)):
                        for j in range(i+1,len(NeighC)):
                            ci = centers[NeighC[i]]
                            cj = centers[NeighC[j]]
                            rij = fc.norm(cj-ci)
                            nij = (cj-ci)/rij
                            if rij <= r0:
                                f_v += 0.1*(rij-r0)*nij
                            else:
                                continue 
        
        F_V.append(f_v)
        
    return np.array(F_V)


def force_vtx_elastic_shear_rest(regions,point_region, ridges, K,A0,G,L, vertices,centers, h,boundary_tissue, bound_down):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (energy gradient descent method)
    
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
    #boundary_tissue = fc.find_boundary_vertices(len(vertices),ridges)
    
    

    for v in range(LV):
        if v not in boundary_tissue:
        
            f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
        
            #Maybe include a regularizing force acting on the cells
        
        
            NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
            for i in range(len(NeighC)):
                for j in range(i+1,len(NeighC)):
                    ci = centers[NeighC[i]]
                    cj = centers[NeighC[j]]
                    rij = fc.norm(cj-ci)
                    nij = (cj-ci)/rij
                    if rij <= r0:
                        f_v += 0.1*(rij-r0)*nij
                    else:
                        continue
        else:
            if v in bound_down:
                f_v = 0.0*np.array([1.,1.]) 
            else:
                f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
                    
                NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
                for i in range(len(NeighC)):
                    for j in range(i+1,len(NeighC)):
                        ci = centers[NeighC[i]]
                        cj = centers[NeighC[j]]
                        rij = fc.norm(cj-ci)
                        nij = (cj-ci)/rij
                        if rij <= r0:
                            f_v += 0.1*(rij-r0)*nij
                        else:
                            continue 
        
        F_V.append(f_v)
        
    return np.array(F_V)



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
    boundary_wound = fc.find_wound_boundary(regions,point_region,wloc)

    for v in range(LV):
        if v not in boundary_tissue:
            if v not in boundary_wound:
                f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
            else:
                f_v = force_vtx_finite_grad_wound(point_region, regions, ridges, vertices, v, K, A0, G, L,Lw,h,wloc, boundary_wound)
            NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
            for i in range(len(NeighC)):
                for j in range(i+1,len(NeighC)):
                    ci = centers[NeighC[i]]
                    cj = centers[NeighC[j]]
                    rij = fc.norm(cj-ci)
                    nij = (cj-ci)/rij
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
    



def smooth_boundary(regions,point_region, ridges, K,A0,G,L, vertices, centers,h):
    LV = len(vertices) #Number of vertices

    F_V = [] # Vector of resulting forces acting on the vertices
    
    
    r0 = np.sqrt(np.mean(A0)/np.pi)
    
    #For vertices, well, this may be a huge mess
    
    #First find the boundary
    boundary_tissue = fc.find_boundary_vertices(len(vertices),ridges)
   

    for v in range(LV):
        if v not in boundary_tissue:
            f_v = 0.0*np.array([1.,1.]) 
        else:
            f_v = force_vtx_finite_grad(point_region, regions, ridges, vertices, v, K, A0, G, L,h,boundary_tissue)
            NeighR, NeighC = fc.find_vertex_neighbour_centers(regions,point_region,v)
            for i in range(len(NeighC)):
                for j in range(i+1,len(NeighC)):
                    ci = centers[NeighC[i]]
                    cj = centers[NeighC[j]]
                    rij = fc.norm(cj-ci)
                    nij = (cj-ci)/rij
                    if rij <= r0:
                        f_v += 0.1*(rij-r0)*nij
                    else:
                        continue
                
        #Maybe include a regularizing force acting on the cells
        F_V.append(f_v)
        
    return np.array(F_V)
    