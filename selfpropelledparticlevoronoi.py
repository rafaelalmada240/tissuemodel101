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
    #Gives all the neighbours of a center in a voronoi tesselation (removes all -1 vertices)    
    point = point_region[center]
    R = regions[point]
    
    for e in R:
        if e==-1:
            R.remove(e)
    return R

def remove_minus(ridges):
    if isinstance(ridges,np.ndarray):
        ridges = ridges.tolist()
    for ridge in ridges:
        for elem in ridge:
            if elem == -1:
                ridges.remove(ridge)
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
    #Calculate the perimeter term - better write it in explicit form to avoid confusion
    for i in range(len(V)):
        P += norm(V[(i+1)%len(V)]-V[i])

    return P

def area_vor(point_region,regions,vertices,i):

    R = find_center_region(regions,point_region,i)
    V = vertices[R]
    A = 0
    #Calculate the area term
    for i in range(len(V)):
        A += norm((np.cross(V[i],V[(i+1)%len(V)])))
        
    return A
        
def nsides_vor(point_region,regions,i):

    R = find_center_region(regions,point_region,i)
    nsides=len(R)
    return nsides

def energy_vor(point_region,regions,ridges, vertices,vertex, K,A0,G,L):

    R,N_c = find_vertex_neighbour_centers(regions, point_region,vertex)
    N_v = find_vertex_neighbour_vertices(ridges,vertex)
    E = 0
    for i in N_c:
        Pi = perimeter_vor(point_region,regions,vertices, i)
        Ai = area_vor(point_region,regions,vertices, i)
        E += K/2*(Ai-A0)**2 + G/2*Pi**2
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
        S = list(set(R_i).intersection(N_v))
        #print('R_i - '+str(R_i)+'\n')
        #print('N_v -'+str(N_v)+'\n')
        
        if len(S)>1:
        
            n_a = -np.cross(z,vertices[S[0]]-vertices[S[1]])[:-1]
        
            v_p0 = vertices[S[0]]-vertices[vertex]
            v_p1 = vertices[S[1]]-vertices[vertex]
            l_0 = norm(v_p0)
            l_1 = norm(v_p1)
            n_p = v_p0/l_0 + v_p1/l_1
        
            Pi = perimeter_vor(point_region,regions,vertices, c_i)
            Ai = area_vor(point_region,regions,vertices, c_i)
            F += K/2*(Ai-A0)*n_a + float(G*Pi+L/2)*n_p
    #for j in N_v:
    #    v = vertices[j]        
    #    edgeV = vertices[vertex] - v
    #    lj = norm(edgeV)
    #    nj = edgeV/lj    
    #    F = L*nj
        
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
            
            rho = 2*r0/(r_nc+1e-1)
            f_c += -5*(rho**2-rho)*l_nc-(r_nc-4*r0)*l_nc
                
        FC.append(f_c)         
        
    return np.array(FC)

     

def force_vtx_elastic(regions,point_region, ridges, K,A0,G,L, vertices,dx):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (energy gradient descent method)
    
    Variables:
    
    regions - (list) set of all the vertices that compose the different regions of the network
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    K, A0, G, L - (float) model parameters 
    vertices - array of coordinate positions for each vertex of the network
    dx - step of gradient descent
    
    Output:

    F_V - array of all forces acting on the vertices

    
    '''
    
    LV = len(vertices) #Number of vertices

    F_V = [] # Vector of resulting forces acting on the vertices
    
    
    #For vertices, well, this may be a huge mess
    new_vertices = vertices
    for v in range(LV):
        
        #f_v = 0.0*np.array([1.,1.])
        Ev = energy_vor(point_region,regions,ridges, vertices,v, K,A0,G,L)
        
        dir = 2*(np.random.rand(2)-0.5)
        ndir = dir/norm(dir)
        
        new_vertices[v] = vertices[v] + ndir*dx
        Ev2 = energy_vor(point_region,regions,ridges, new_vertices,v, K,A0,G,L)
        f_v = -(Ev2-Ev)*ndir/(dx)
        
        new_vertices = vertices
        F_V.append(f_v)
        
    return np.array(F_V)


def force_vtx_elastic_v1(regions,point_region, ridges, K,A0,G,L, vertices):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (finite network difference approach)
    
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
    
    
    #For vertices, well, this may be a huge mess
    for v in range(LV):
        
        f_v = 0.0*np.array([1.,1.])
        Ev = energy_vor(point_region,regions,ridges, vertices,v, K,A0,G,L)
        NeighC_V, NeighV_V = find_vertex_neighbour(regions,point_region,ridges,v)
        
        for v2 in NeighV_V:
            
            r2 = vertices[v2]
                    
            #Direction
            edgeV = r2 - vertices[v]
            lv = norm(edgeV)
            nv = edgeV/(lv+0.1)
                    
            #Calculate the local energy for the neighbouring vertices
            Ev2 = energy_vor(point_region,regions,ridges, vertices,v2, K,A0,G,L)
                
            #Implement the force as an energy gradient
            f_v += (Ev-Ev2)*nv/(lv+0.01)
            
        F_V.append(f_v)
        
    return np.array(F_V)

def cells_avg_vtx(regions,point_region,cells,vertices):
    for i in range(len(cells)):
        Neigh_c = find_center_region(regions,point_region, i)
        avg_vc = np.mean(vertices[Neigh_c],0)
        cells[i] = avg_vc
        
    return cells
        
    

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
        
        f_v = force_vtx(point_region,regions,ridges,vertices,v,K,A0,G,L)#0.0*np.array([1.,1.])    
            
        F_V.append(f_v)
        
    return np.array(F_V)
            
    
## T1 Transitions

def T1_rotations(vertices, v_orig,v_min):
    '''Rotates the coordinates of an edge undergoing a T1 transition
    Input:
    vertices - set of coordinate values for each vertex
    (v_orig,v_min) - edge in the iteration that is undergoing a T1 transition
    
    Output:
    (r_vorig, r_vmin) - new coordinate values of the rotated edge
    '''
    th = np.pi/2
    
    r_vorig = vertices[v_orig]
    r_vmin = vertices[v_min]
    R_p = np.array([[0,-1],[1,0]])
    R_m = np.array([[0,1],[-1,0]])
                    
    diff_v_vm = r_vorig-r_vmin
    new_diff_v_vm_coord_1 = r_vmin + np.sqrt(3)/2*(R_p.dot(diff_v_vm)+1/2*diff_v_vm)
    new_diff_v_vm_coord_2 = r_vmin + np.sqrt(3)/2*(R_m.dot(diff_v_vm)+1/2*diff_v_vm)
            
        
    # Assign each vertex to each of the ends of the rotated vector
        
    p_1 = np.random.rand()
        
    if p_1 > 0.5:
        r_vorig = new_diff_v_vm_coord_1
        r_vmin = new_diff_v_vm_coord_2
    else:
        r_vorig = new_diff_v_vm_coord_2
        r_vmin = new_diff_v_vm_coord_1
        
    return r_vorig, r_vmin

def T1_len_neigh_after_rot(vertices, v_orig,neigh_min,v_min,neigh_orig):
    '''
    During a T1 transition finds the vertices for the new neighbourhood for the edge undergoing a T1 transition
    Input:
    vertices - set of coordinate values for each vertex
    (v_orig,v_min) - edge in the iteration that is undergoing a T1 transition
    neigh_min - neighbourhood of v_min
    neigh_orig - neighbourhood of v_orig
    
    Output:
    new_neighv,new_neighm - Neighbourhoods of vorig and vmin that are not vorig or vmin, respectively
    lvv_min, lvm_min - the index of the vertices in new_neighv and new_neighm that will be swapped between the two neighbourhoods, respectively
    '''
    
    
    list_len_vm = []
    new_neighm = []
                        
    for v in neigh_min:
        if v != v_orig:
            #For this one we don't need to remove previous vertices
            deltaxorig = vertices[v_orig]-vertices[int(v)]
            list_len_vm.append(norm(deltaxorig))
            new_neighm.append(int(v))
        
    list_len_v = []
    new_neighv = []
    for v in neigh_orig:
        if v != v_min:
            deltaxmin = vertices[v_min]-vertices[v]
            list_len_v.append(norm(deltaxmin))
            new_neighv.append(v)
            
            lvm_min = np.argmin(np.array(list_len_vm))
            #lvv_min = np.argmin(np.array(list_len_v))
                        
    if lvm_min == 0:
        lvv_min = 1
    else: 
        lvv_min = 0
    return new_neighv,new_neighm, lvv_min, lvm_min

def T1_change_edges(ridges,vertices_neigh, vertices_neigh_min,i,v_min):
    loc_ridges = np.where(ridges == i)[0]
            
    loc_neigh_not_vm = np.where(np.array(vertices_neigh)!=v_min)[0]

    skip_parameter = int(0)
    for j in range(len(loc_ridges)):
                
        if v_min in ridges[loc_ridges[j]]:
            skip_parameter += int(1)
            continue
        else:
            js = int(j-skip_parameter)
            ridges[loc_ridges[j]]= [vertices_neigh[loc_neigh_not_vm[js]],i]                    
    loc_ridges = np.where(ridges == v_min)[0]
            
    loc_neigh_not_i = np.where(np.array(vertices_neigh_min)!= i)[0]
            
    skip_parameter = int(0)
                        
    for j in range(len(loc_ridges)):
        if i in ridges[loc_ridges[j]]:
            skip_parameter+=int(1)
            continue
        else:
                                
            js = int(j-skip_parameter)
            ridges[loc_ridges[j]]= [vertices_neigh_min[loc_neigh_not_i[js]],v_min]
                                
    return ridges
    

def T1transition(vertices, ridges, regions, point_region,thresh_len):
    
    ''' This function runs through all the interior vertices on the network and does topological rearrangements (T1 transitions) 
    if the length of the edges are below a certain threshold. The overall effect on the network should be similar to an viscoelastic relaxation 
    
    Variables:
    
    vertices - list of coordinate positions for each vertex of the network
    ridges - set of all edges in the network
    regions - set of all the vertices that compose the different regions of the network
    point_region - set of the different regions of the network
    thresh_len - parameter that determines the transition length
    
    Output
    
    Updated versions of vertices, ridges and regions
    transition_counter = number of T1 transitions at this iteration of the implementation
    
    '''
    
    #make a list of vertices to run through in the next cycle
    
    transition_counter = 0
    
    vertex_list = []
    
    for k in range(len(vertices)):
        vertex_list.append(k)
        
    
    for i in vertex_list:
        
        #First check if vertex is not empty or blown up, otherwise skip
        if np.isnan(np.sum(vertices[i])):
            continue
        else:
        
            #Find all neighbouring vertices of vertex i
        
            vertices_neigh = find_vertex_neighbour_vertices(ridges,i)
        
            # list of neighbouring vertices lengths
            list_len = []
            list_neigh_v_not_excluded = [] #List of the neighbours of v that have not been excluded yet
            for v in vertices_neigh:
                
                #Only include the vertices that have not been in cycle yet
                if v in vertex_list:
                    list_neigh_v_not_excluded.append(int(v))
                    deltax = vertices[i]-vertices[int(v)]
                    list_len.append(np.sqrt(np.sum(deltax**2)))
        
            if len(list_neigh_v_not_excluded) > 2:
            
                # Find closest neighbouring vertices
        
                loc_v_min = np.argmin(list_len)
                lv_min = list_len[loc_v_min]
                v_min = int(list_neigh_v_not_excluded[loc_v_min])
        
                # Find neighbours of closest vertex to vertex i
        
                vertices_neigh_min = find_vertex_neighbour_vertices(ridges,v_min)
                
                #Only do this for vertices with 3 neighbours
            
                if len(vertices_neigh_min) > 2:
                    #Do rotation of the vertex edge along 90 degrees (The geometric part is okay, not likely to change much)
                    if lv_min < thresh_len: 
                        vertices[i], vertices[v_min] = T1_rotations(vertices,i,v_min)
                        
                        #list of neighbouring vertice lengths in the closest neighbouring vertex for rotated vector only assign vertices in vertex list        
                        new_neighv, new_neighvm, lvv_min, lvm_min = T1_len_neigh_after_rot(vertices,i,vertices_neigh_min,v_min,list_neigh_v_not_excluded)        
                        vertices_neigh_min = list(set(vertices_neigh_min).difference({new_neighvm[lvm_min]}).union({new_neighv[lvv_min]}))
                        vertices_neigh = list(set(vertices_neigh).difference({new_neighv[lvv_min]}).union({new_neighvm[lvm_min]}))
                                    
                        #Wait, you actually need to change the ridges as well, hopefully this works                        
                        ridges = T1_change_edges(ridges,vertices_neigh,vertices_neigh_min,i,v_min)
                                   
                        #Oh god, I forgot about the connections between the vertices and the centers which also change                       
                        #For vertex i
                        regions_neigh_v = []
                        for v in vertices_neigh:
                            ff = find_vertex_neighbour_centers(regions,point_region,v)
                            if ff is None:
                                continue
                            else:
                                region, center = ff
                                regions_neigh_v.append(region)
                        region_v = []
                        
                        for k in range(len(regions_neigh_v)):
                            region_v.append(list(set(regions_neigh_v[k]).intersection(regions_neigh_v[(k+1)%len(regions_neigh_v)])))
                            
                        region_v_true = []
                        
                        for r in region_v:
                            region_v_true = list(set(region_v_true).union(r))
                            
                            
                        for r in region_v_true:
                            if i in regions[r]:
                                continue
                            else:
                                
                                if v_min in regions[r]:
                                    regions[r].append(i)
                                    regions[r].remove(v_min)
                                
                                
                        #For neighbouring vertex v_min
                        regions_neigh_vmin = []
                        for v in vertices_neigh_min:
                            ff = find_vertex_neighbour_centers(regions,point_region,v)
                            if ff is None:
                                continue
                            else:
                                region, center = ff
                                regions_neigh_vmin.append(region)
                            
                        region_vmin = []
                        
                        for k in range(len(regions_neigh_vmin)):
                            region_vmin.append(list(set(regions_neigh_vmin[k]).intersection(regions_neigh_vmin[(k+1)%len(regions_neigh_vmin)])))
                            
                        region_vmin_true = []
                        
                        for r in region_vmin:
                            region_vmin_true = list(set(region_vmin_true).union(r))
                            
                        for r in region_vmin_true:
                            if v_min in regions[r]:
                                continue
                            else:
                                if i in regions[r]:
                                    regions[r].remove(i)
                                    regions[r].append(v_min)
                                    
                        
                        transition_counter += 1
                        #with open('Transitions'+str(i)+'.txt','a') as text:
                            #text.write(' Neighs of v prior to transition: ')
                            #text.write(str(list_neigh_v_not_excluded)+ '\n')
                            #text.write('V_min: '+str(v_min)+'\n') 
                            #text.write(' Neighs of v_min prior to transition: ')   
                            #text.write(str(vertices_neigh_min_o)+'\n') 
                            #text.write('LENGTHS OF NEIGH_v TO V_MIN: '+str(list_len_v)+'\n')
                            #text.write('LENGTHS OF NEIGH_vmin TO V: '+str(list_len_vm)+'\n')
                        
                        
                            
                        
                            
                                
                    vertex_list.remove(i)
                    
            elif len(list_neigh_v_not_excluded) <= 2:
                continue

    return ridges, vertices, regions, transition_counter


## Depricated functions (or not necessary anymore)

def path_region(regions, point_region, vertices, center):
    #Seems to not be necessary
    R = find_center_region(regions,point_region, center)
    dist_mat = np.zeros((len(R),len(R)))
    for i in range(len(R)):
        for j in range(len(R)):
            dist_mat[i,j] = np.sum((vertices[R[i]]-vertices[R[j]])**2)**0.5
            
    R_new = []
    R_new.append(R[0])
    
    for i in range(len(R)-1):
        j = np.where(np.array(R) == R_new[i])[0][0]
        nl = []
        for k in range(len(R)):
            if R[k] not in R_new:
                nl.append(k)
        v = np.argmin(dist_mat[j,nl])
        R_new.append(R[v])
        
    return R

def force_vor_elastic(regions,point_region, ridges,boundary, Kp,r0, cells, vertices):
    '''
    
    Calculates the force in all of the vertices according to a viscoelastic network model
    
    Variables:
    
    regions - (list) set of all the vertices that compose the different regions of the network
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    boundary - (list) set of all boundary vertices
    Kp - (float) elastic parameter 
    r0 - (float) natural spring length
    cells - array of coordinates for all cell centers of the network
    vertices - array of coordinate positions for each vertex of the network
    
    Output:
    F_C - array of all forces acting on the cell centres
    F_V - array of all forces acting on the vertices
    mean_dist_v - (float) average distance between vertices
    
    '''
    
    LC = len(cells)
    LV = len(vertices)
    F_C = [] # Vector of resulting forces acting on the centers
    F_V = [] # Vector of resulting forces acting on the vertices
    
    dist_v_list = []
    
    
    for c in range(LC):
        Neigh_c = find_center_region(regions,point_region, c)
        
        f_c = 0.0*np.array([1.,1.])
        for n_c in Neigh_c:
            v_c = vertices[n_c]
            avg_vc = 0
            if any(v_c < 5) and any(v_c > -5):
                r_vec_center = cells[c] - v_c
                avg_vc += v_c/len(Neigh_c)
                abs_rvc = norm(r_vec_center)
                n_rvc = r_vec_center/abs_rvc
                f_c += -Kp*(abs_rvc-r0)*n_rvc
            #r_dev = cells[c] - avg_vc
            #abs_dev = norm(r_dev)
            #n_dev = r_dev/abs_dev
            #f_c += Kp*(abs_dev)*n_dev
        F_C.append(f_c)
    
    for v in range(LV):
        f_v = 0.0*np.array([1.,1.])
        if True:#v in boundary:
            #f_v = 0.0*np.array([1.,1.])
        #else:
             
            NeighC_V, NeighV_V = find_vertex_neighbour(regions,point_region,ridges,v)
            
            #Vertices in the boundary can have 2 neigbours or less, so fixing them is equivalent to fixing the boundary conditions in a certain way.
            
            for nc_v in NeighC_V:
                c_v = cells[nc_v]
                r_center_vertex = vertices[v] - c_v
                abs_rcv = norm(r_center_vertex)
                n_rcv = r_center_vertex/abs_rcv
            
                f_v += -Kp*(abs_rcv - r0)*n_rcv
        
            for nv_v in NeighV_V:
                v_v = vertices[nv_v]
                if any(v_v < 5) and any(v_v > -5):
                    edgeV = vertices[v] - v_v
                    abs_vv = norm(edgeV)
                    dist_v_list.append(abs_vv)
                    n_vv = edgeV/abs_vv
            
                    f_v += -Kp*(abs_vv - r0)*n_vv
            
        F_V.append(f_v)
        
    return np.array(F_C), np.array(F_V), np.mean(dist_v_list)

def phi_n(l_n):
    return l_n*np.tan(np.pi/l_n)

def find_boundary_vertices(n_vertices,ridges):
    '''
    This function finds all vertices in the boundary of the network, under the initial assumption that vertices have 3 connections in planar graphs if they are in a bulk and 2 edges or less if they are in the boundary
    
    Variables:
    
    n_vertices - is the number of vertices in the network (the cardinality of the vertex set)
    
    ridges - is the set of all edges in the network
    
    Output:
    
    Bound_set - the set containing all vertices that are in the boundary, this list doesn't contain all vertices in the boundary, but it ensures all vertices in the list are in the boundary.
    '''
    vertex_list = []
    
    for k in range(n_vertices):
        vertex_list.append(k)
    Bound_set = []
    Bound_set1 = []
    Bound_set_neighbours = []
    for v in vertex_list:
        Neigh_V = find_vertex_neighbour_vertices(ridges,v)
        if len(Neigh_V) < 3:
            Bound_set.append(v)
            Bound_set1.append(v)
            Bound_set_neighbours.append(Neigh_V)
        vertex_list.remove(v)
        
    for i in range(len(Bound_set1)):
    #    for j in range(i+1,len(Bound_set1)):
        neigh1 = Bound_set_neighbours[i]
    #        neigh2 = Bound_set_neighbours[j]
    #        vertex_new = list(set(neigh1).intersection(neigh2))
        for b in neigh1:
            #vertex_new:
            if b not in Bound_set:
                Bound_set.append(b)
                    
    return Bound_set
                
# Probably not necessary
def remove_empty(regions):
    for region in regions:
        if len(region)==0:
            regions.remove(region)
    return regions