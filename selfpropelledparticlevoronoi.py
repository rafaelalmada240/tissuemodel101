import numpy as np

'''
    regions - (list) set of all the vertices that compose the different regions of the network
    point_region - (list) set of the different regions of the network
    ridges - (list) set of all edges in the network
    cells - array of coordinates for all cell centers of the network
    vertices - array of coordinate positions for each vertex of the network

'''

def norm(vec):
    return np.sqrt(np.sum(vec**2))

def path_region(regions, point_region, vertices, center):
    R = find_center_region(regions,point_region, center)
    dist_mat = np.zeros((len(R),len(R)))
    for i in range(len(R)):
        for j in range(len(R)):
            dist_mat[i,j] = np.sum((vertices[R[i]]-vertices[R[j]])**2)**0.5
            
    R_new = []
    R_new.append(R[0])
    
    for i in range(len(R)-1):
        j = np.where(np.arrayR == R_new[i])[0][0]
        nl = []
        for k in range(len(R)):
            if R[k] not in R_new:
                nl.append(k)
        v = np.argmin(dist_mat[j,nl])
        R_new.append(R[v])
        
    return R
        

#Functions to find locations of vertices and centers

def find_center_region(regions,point_region,center):
    #Gives all the neighbours of a center in a voronoi tesselation (removes all -1 vertices)
    R = regions[point_region[center]]
    for e in R:
        if e==-1:
            R.remove(e)
    return R


def remove_minus(ridges):
    ridges = list(ridges)
    for ridge in ridges:
        for elem in ridge:
            if elem == -1:
                ridges.remove(ridge)
    return np.array(ridges)

def find_vertex_neighbour_vertices(ridges,vertex):
    
    #Gives all the vertex neighbours of a vertex in a voronoi tesselation (removes all -1 vertices)

    list_vertex_neigh = []
        
    for ridge in ridges:
        if vertex in list(ridge):  
            for elem in list(ridge):
                if (elem != vertex) and (elem != -1):
                    list_vertex_neigh.append(int(elem))
    
    return np.sort(list_vertex_neigh)

def find_vertex_neigbour_centers(regions, point_region,vertex):
    list_regions = []
    list_centers = []
    i = 0
    for region in regions:
        if len(region)>0:
            if vertex in region:
                list_regions.append(i)
        i+=1       
    for region in list_regions:
        loc_points = np.where(point_region==region)
        list_centers.append(loc_points[0][0])
        
        
        
        return list_regions, np.array(list_centers)
    
def find_vertex_neighbour(regions, point_region, ridges,vertex):
    
    #Gives all the neighbours of a vertex in a voronoi tesselation (removes all -1 vertices)
    #print(regions)
    #print(point_region)
    #print(vertex)
    #print(find_vertex_neigbour_centers(regions, point_region,vertex))
    
    ff = find_vertex_neigbour_centers(regions, point_region,vertex)
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
        A, L_c = find_vertex_neigbour_centers(regions,point_region,v)
        List_centers = list(set(List_centers).union(L_c))
    return List_centers

## Functions to calculate perimeter and area terms for a given vertex

def perimeter_vor(point_region,regions,vertices,i):
    R = path_region(point_region,regions,vertices,i)
    V = vertices[R]
    P = np.sum(np.sqrt(np.sum(np.diff(V,axis=0)**2,1)))+np.sqrt(np.sum((V[-1]-V[0])**2))
    #Calculate the perimeter term
    return P

def area_vor(point_region,regions,vertices,i):
    R = path_regions(point_region,regions,vertices,i)
    V = vertices[R]
    A = 0
    #Calculate the area term
    for i in range(len(V)):
        A += np.cross(np.cross(v[i],v[i+1]))
    return A
        
def nsides_vor(point_region,region,i):
    R = find_center_region(point_region,region,i)
    nsides=len(R)
    return nsides

def energy_vor(P,P0,Kp):
    
    E = np.sum(Kp*(P-P0)**2)
    return E


## Different force models 




def thermalize(regions, point_region,cells,vel,r0):
    LC = len(cells)
    #print(LC)
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
            r_dev = cells[c] - avg_vc
            abs_dev = norm(r_dev)
            n_dev = r_dev/abs_dev
            f_c += Kp*(abs_dev)*n_dev
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

def force_vtx_elastic(regions,point_region, ridges,boundary, Kp,r0, cells, vertices):
    
    '''
    
    Calculates the force in all of the vertices according to the vertex model (not yet implemented)
    
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
            if any(v_c < 5) and any(v_c > -5):
                r_vec_center = cells[c] - v_c
                abs_rvc = norm(r_vec_center)
                n_rvc = r_vec_center/abs_rvc
                f_c += -Kp*(abs_rvc-r0)*n_rvc
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
                
            
    
## T1 Transitions

def T1transition(vertices, ridges, regions, point_region,thresh_len):
    
    ''' This function runs through all the interior vertices on the network and does topological rearrangements (T1 transitions) if the length of the edges are below a certain threshold. The overall effect on the network should be similar to an viscoelastic relaxation 
    
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
        
        if np.isnan(np.sum(vertices[i])):
            continue
        else:
        
        
            #Find all neighbouring vertices of vertex i
        
            vertices_neigh = find_vertex_neighbour_vertices(ridges,i)
        
            # list of neighbouring vertices lengths
            list_len = []
            list_neigh_v_not_excluded = []
            for v in vertices_neigh:
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
            
                if len(vertices_neigh_min) > 2:
        
                    #Do rotation of the vertex edge along 90 degrees
                    th = np.pi/2
                    R_p = np.array([[0,-1],[1,0]])
                    R_m = np.array([[0,1],[-1,0]])
        
        
                    new_diff_v_vm_coord_1 = 0
                    new_diff_v_vm_coord_2 = 0
                    
        
                    if lv_min < thresh_len:
                        diff_v_vm = vertices[i]-vertices[v_min]
                        new_diff_v_vm_coord_1 = vertices[v_min] + np.sqrt(3)/2*(R_p.dot(diff_v_vm)+1/2*diff_v_vm)
                        new_diff_v_vm_coord_2 = vertices[v_min] + np.sqrt(3)/2*(R_m.dot(diff_v_vm)+1/2*diff_v_vm)
                        
                        #with open('vertextransitions.txt','a') as text:
                    
                        #    text.write('Old length: ')
                    
                        #    text.write(str(lv_min)+' and the ')
                    
                        #    text.write('threshold: ')
                    
                        #    text.write(str(thresh_len)+'\n')
                            
                        #    transition_counter += 1
            
        
                        # Assign each vertex to each of the ends of the rotated vector
        
                        p_1 = np.random.rand()
        
                        if p_1 > 0.5:
                            vertices[i] = new_diff_v_vm_coord_1
                            vertices[v_min] = new_diff_v_vm_coord_2
                        else:
                            vertices[i] = new_diff_v_vm_coord_2
                            vertices[v_min] = new_diff_v_vm_coord_1
            
        
                        #list of neighbouring vertice lengths in the closest neighbouring vertex for rotated vector only assign vertices in vertex list
        
                        list_len_vm = []
                        new_neighvm = []
                        
                        for v in vertices_neigh_min:
                            if v != i:
                                #For this one we don't need to remove previous vertices
                                deltax = vertices[i]-vertices[int(v)]
                                list_len_vm.append(np.sqrt(np.sum(deltax**2)))
                                new_neighvm.append(int(v))
        
                        list_len_v = []
                        new_neighv = []
                        for v in list_neigh_v_not_excluded:
                            if v != v_min:
                                deltax = vertices[v_min]-vertices[v]
                                list_len_v.append(np.sqrt(np.sum(deltax**2)))
                                new_neighv.append(v)
        
                        lvm_min = np.argmin(np.array(list_len_vm))
                        lvv_min = np.argmin(np.array(list_len_v))
                
        
        
                        vertices_neigh_min = list(vertices_neigh_min)
                        vertices_neigh = list(vertices_neigh)
        
                        vertices_neigh_min.append(new_neighv[lvv_min])
                        vertices_neigh.append(new_neighvm[lvm_min])
        
                        vertices_neigh.remove(new_neighv[lvv_min])
                        vertices_neigh_min.remove(new_neighvm[lvm_min])
            
                        #Wait, you actually need to change the ridges as well, hopefully this works
            
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
                                
                        
                        #Oh god, I forgot about the connections between the vertices and the centers which also change
                        
                        #For vertex i
                        regions_neigh_v = []
                        for v in vertices_neigh:
                            ff = find_vertex_neigbour_centers(regions,point_region,v)
                            if ff is None:
                                continue
                            else:
                                region, center = ff
                                regions_neigh_v.append(region)
                            
                        #print(regions_neigh_v)
                        region_v = []
                        
                        for k in range(len(regions_neigh_v)):
                            region_v.append(list(set(regions_neigh_v[k]).intersection(regions_neigh_v[(k+1)%len(regions_neigh_v)])))
                            
                        region_v_true = []
                        
                        for r in region_v:
                            region_v_true = list(set(region_v_true).union(r))
                            
                        #print(region_v_true)
                            
                        for r in region_v_true:
                            if i in regions[r]:
                                continue
                            else:
                                #print(i)
                                #print(v_min)
                                #print(regions[r])
                                
                                if v_min in regions[r]:
                                    regions[r].append(i)
                                    regions[r].remove(v_min)
                                
                                #print(regions[r])
                                
                        #For neighbouring vertex v_min
                        regions_neigh_vmin = []
                        for v in vertices_neigh_min:
                            ff = find_vertex_neigbour_centers(regions,point_region,v)
                            if ff is None:
                                continue
                            else:
                                region, center = ff
                                regions_neigh_vmin.append(region)
                            
                        region_vmin = []
                        
                        for k in range(len(regions_neigh_vmin)):
                            region_vmin.append(list(set(regions_neigh_vmin[k]).intersection(regions_neigh_vmin[(k+1)%len(regions_neigh_v)])))
                            
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
                                
                        
                        
                            
                        
                            
                                
                    vertex_list.remove(i)
            elif len(list_neigh_v_not_excluded) <= 2:
                continue

    return ridges, vertices, regions, transition_counter