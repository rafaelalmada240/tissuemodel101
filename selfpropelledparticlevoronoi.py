import numpy as np

'''
    regions - (list) set of all the vertices that compose the different regions of the network
    Some regions are empty, but it is best not to remove them for consistency
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    cells - array of coordinates for all cell centers of the network
    vertices - array of coordinate positions for each vertex of the network

'''

def norm(vec):
    return np.sqrt(np.sum(vec**2))

#Auxiliary functions to find locations of vertices and centers

def find_center_region(regions,point_region,center):
    '''Gives all the neighbours of a center in a voronoi tesselation (removes all -1 vertices)'''
    point = point_region[center]
    R = regions[point]
    
    for e in R:
        if e==-1:
            R.remove(e)
    return R

def remove_minus(ridges):
    if isinstance(ridges,np.ndarray):
        ridges = ridges.tolist()
        
    index_to_remove = []
    for ridge in ridges:
        for elem in ridge:
            if elem == -1:
                index_to_remove.append(ridge)
    
    for j in index_to_remove:
        ridges.remove(j)
    return ridges


def find_vertex_neighbour_vertices(ridges,vertex):
    
    '''
    The function receives as input the list of ridges (equivalent to the set of edges E) and 
    an integer corresponding to the vertex v_i, and return a list of vertices v_j, such that (v_i, v_j) is in E
    '''

    list_vertex_neigh = []
        
    for ridge in ridges:
        if vertex in list(ridge):  
            for elem in list(ridge):
                if (elem != vertex) and (elem != -1):#This condition ensures that we don't include either the vertex v_i or the vertex at infinity
                    list_vertex_neigh.append(int(elem))
    
    return np.sort(list_vertex_neigh)

def find_vertex_neighbour_centers(regions, point_region,vertex):
    
    '''
    Inputs: 
    regions - (list) set of all the vertices that compose the different regions of the network
    Some regions are empty, but it is best not to remove them for consistency
    point_region - (list) set of the different regions of the network
    vertex - (int) a specific vertex
    
    Outputs:
    list_regions, list_centers - list of regions and corresponding centers that are neighbouring a vertex.
    
    '''
    list_regions = []
    list_centers = []
    i = 0

        
    for i in range(len(regions)):
    #Only consider neighbouring regions that form polygons
        if vertex in regions[i]:
            list_regions.append(i) 
            loc_points = np.where(np.array(point_region)==i)
            list_centers.append(loc_points[0][0])
             
    return list_regions, np.array(list_centers)

def find_vertex_neighbour(regions, point_region, ridges,vertex):
    
    '''
    Inputs: 
    regions - (list) set of all the vertices that compose the different regions of the network
    Some regions are empty, but it is best not to remove them for consistency
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    vertex - (int) a specific vertex
    
    Outputs:
    list_vertex neigh, list_centers - list of centers and vertices that are neighbouring a vertex.
    
    '''
    
    #Gives all the neighbours of a vertex in a voronoi tesselation (removes all -1 vertices)
    
    ff = find_vertex_neighbour_centers(regions, point_region,vertex)
    if ff is None:
        list_centers = []
    else:
        list_regions, list_centers = ff 
    list_vertex_neigh = find_vertex_neighbour_vertices(ridges,vertex)
    
    
    
    return np.sort(list_centers), np.sort(list_vertex_neigh)    

def find_center_neighbour_center(regions,point_region,center):
    
    # Find all neighbouring cells for a given cell i
    List_centers = []

    R = find_center_region(regions, point_region, center)
    for v in R:
        A, L_c = find_vertex_neighbour_centers(regions,point_region,v)
        List_centers = list(set(List_centers).union(L_c))
    List_centers.remove(center)
    return List_centers



## Functions to calculate perimeter and area terms for a given vertex

def perimeter_vor(point_region,regions,vertices,i):

    R = find_center_region(regions,point_region,i)
    V = vertices[R]
    P = 0
    #Calculate the perimeter term
    if len(R)>2:
        for i in range(len(V)):
            P += norm(V[(i+1)%len(V)]-V[i])
    return P

def area_vor(point_region,regions,vertices,i):

    R = find_center_region(regions,point_region,i)
    V = vertices[R]
    A = 0
    #Calculate the area term
    if len(R)>2:
        for i in range(len(V)):
            A += 1/2*norm((np.cross(V[i],V[(i+1)%len(V)])))        
    return A
        
def nsides_vor(point_region,regions,i):

    R = find_center_region(regions,point_region,i)
    nsides=len(R)
    return nsides

def energy_vor(point_region,regions,ridges, vertices,vertex, K,A0,G,L):

    R,N_c = find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = find_vertex_neighbour_vertices(ridges,vertex)
    E = 0
    
    #Sum over all cells containing the vertex
    for i in N_c:
        Pi = perimeter_vor(point_region,regions,vertices, i)
        Ai = area_vor(point_region,regions,vertices, i)
        E += K/2*(Ai-A0)**2 + G/2*Pi**2
        
    #Sum over all adjacent vertices
    for j in N_v:
        v = vertices[j]        
        edgeV = vertices[vertex] - v
        lj = norm(edgeV)    
        E += L*lj
        
    return E


def force_vtx(point_region, regions, ridges, vertices, vertex, K, A0, G, L):
    
    R,N_c = find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = find_vertex_neighbour_vertices(ridges,vertex)
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
                l_0 = norm(v_p0)
                l_1 = norm(v_p1)
                n_p = v_p0/l_0 - v_p1/l_1
        
                Pi = perimeter_vor(point_region,regions,vertices, c_i)
                Ai = area_vor(point_region,regions,vertices, c_i)
                
                print(Pi)
                print(Ai)
                F += K/2*(Ai-A0)*n_a - float(G*Pi+L/2)*n_p
        
    return F
## Different force models 

def thermalize(regions, point_region,cells,vel,r0):
    LC = len(cells)
    FC = []
    
    for c in range(LC):
        xy_c = cells[c]
        neigh_c = find_center_neighbour_center(regions,point_region,c)
        f_c = 0.0*np.array([1.,1.])#+(np.random.rand(2)-0.5)*2
        for n_c in neigh_c:
            xy_n = cells[n_c]
            v_nc = xy_c-xy_n
            r_nc = norm(v_nc)
            l_nc =v_nc/(r_nc+1e-1)
            
            #Lennard Jones
            
            rho = r0/(r_nc+1e-1)
            f_c += -5*(rho**2-rho)*l_nc-(r_nc-2*r0)*l_nc
                
        FC.append(f_c)         
        
    return np.array(FC)

     

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
    
    new_vertices = np.array(vertices)
    old_vertices = np.array(vertices)
    
    r0 = np.sqrt(A0/np.pi)
    
    #For vertices, well, this may be a huge mess

    for v in range(LV):
        new_vertices1x = np.array(new_vertices)
        new_vertices2x = np.array(new_vertices)
        
        new_vertices1y = np.array(new_vertices)
        new_vertices2y = np.array(new_vertices)
        
        f_v = 0.0*np.array([1.,1.])
        Ev = energy_vor(point_region,regions,ridges, old_vertices,v, K,A0,G,L)
        
        n1 = np.array([1,0])
        n2 = np.array([0,1])
        
        new_vertices1x[v] = new_vertices1x[v] - h*n1
        new_vertices2x[v] = new_vertices2x[v] + h*n1
        
        new_vertices1y[v] = new_vertices1y[v] - h*n2
        new_vertices2y[v] = new_vertices2y[v] + h*n2
        
        Ev1x = energy_vor(point_region,regions,ridges, new_vertices1x,v, K,A0,G,L)
        Ev2x = energy_vor(point_region,regions,ridges, new_vertices2x,v, K,A0,G,L)
        
        Ev1y = energy_vor(point_region,regions,ridges, new_vertices1y,v, K,A0,G,L)
        Ev2y = energy_vor(point_region,regions,ridges, new_vertices2y,v, K,A0,G,L)
        
        dEdx = 0.5*(Ev2x-Ev1x)/h
        dEdy = 0.5*(Ev2y-Ev1y)/h
        
        f_v = -(dEdx*n1 + dEdy*n2)
        
        #Maybe include a regularizing force acting on the cells
        
        new_vertices = np.array(old_vertices)
        
        NeighR, NeighC = find_vertex_neighbour_centers(regions,point_region,v)
        for i in range(len(NeighC)):
            for j in range(i+1,len(NeighC)):
                ci = centers[NeighC[i]]
                cj = centers[NeighC[j]]
                rij = norm(cj-ci)
                nij = (cj-ci)/rij
                if rij <= r0:
                    f_v += (rij-r0)*nij
                else:
                    continue
        
        F_V.append(f_v)
        
    return np.array(F_V)


def cells_avg_vtx(regions,point_region,cells,vertices):
    for i in range(len(cells)):
        Neigh_c = find_center_region(regions,point_region, i)
        avg_vc = np.mean(vertices[Neigh_c],0)
        cells[i] = avg_vc
        
    return cells
            
            
#Depricated, no longer in use
def force_vtx_elastic_v2(regions,point_region, ridges, K,A0,G,L, vertices):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (analytical expression)
    
    Variables:
    
    regions - (list) set of all the vertices that compose the different regions of the network
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    K, A0, G, L - (float) model parameters 
    vertices - array of coordinate positions for each vertex of the network
    
    Output:
    F_V - array of all forces acting on the vertices
    
    '''
    LV = len(vertices) #Number of vertices
    F_V = [] # Vector of resulting forces acting on the vertices
    
    dist_v_list = []
    
    #For vertices, well, this may be a huge mess
    for v in range(LV):
        
        f_v = force_vtx(point_region,regions,ridges,vertices,v,K,A0,G,L)   
            
        F_V.append(f_v)
        
    return np.array(F_V)
    



