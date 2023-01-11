import numpy as np

def find_center_region(vor_regions,vor_point_region,center):
    #Gives all the neighbours of a center in a voronoi tesselation (removes all -1 vertices)
    R = vor_regions[vor_point_region[center]]
    for e in R:
        if e==-1:
            R.remove(e)
    return R


def remove_minus(vor_ridge):
    vor_ridge = list(vor_ridge)
    for ridge in vor_ridge:
        for elem in ridge:
            if elem == -1:
                vor_ridge.remove(ridge)
    return np.array(vor_ridge)

def find_vertex_neighbour_vertices(vor_ridges,vertex):
    
    #Gives all the vertex neighbours of a vertex in a voronoi tesselation (removes all -1 vertices)

    list_vertex_neigh = []
        
    for ridge in vor_ridges:
        if vertex in ridge:            
            for elem in ridge:
                if (elem != vertex) and (elem != -1):
                    list_vertex_neigh.append(int(elem))
    
    return np.sort(list_vertex_neigh)

def find_vertex_neigbour_centers(vor_regions, vor_point_region,vertex):
    list_regions = []
    list_centers = []
    i = 0
    for region in vor_regions:
        if vertex in region:
            list_regions.append(i)
        i+=1       
    for region in list_regions:
        loc_points = np.where(vor_point_region==region)
        list_centers.append(loc_points[0][0])
        
        return list_regions, np.array(list_centers)
    
def find_vertex_neighbour(vor_regions, vor_point_region, vor_ridges,vertex):
    
    #Gives all the neighbours of a vertex in a voronoi tesselation (removes all -1 vertices)
    
    list_regions, list_centers = find_vertex_neigbour_centers(vor_regions, vor_point_region,vertex)
    list_vertex_neigh = find_vertex_neighbour_vertices(vor_ridges,vertex)
    
    return np.sort(list_centers), np.sort(list_vertex_neigh)    


def perimeter_vor(vor,i):
    k = vor.point_region[i]
    R = vor.regions[k]
    for e in R:
        if e==-1:
            R.remove(e)
    V = vor.vertices[R]
    
    P = np.sum(np.sqrt(np.sum(np.diff(V,axis=0)**2,1)))+np.sqrt(np.sum((V[-1]-V[0])**2))
    
    return P

def nsides_vor(vor,i):
    k = vor.point_region[i]
    R = vor.regions[k]
    for e in R:
        if e==-1:
            R.remove(e)
    nsides=len(R)
    
    return nsides

def energy_vor(P,P0,Kp):
    
    E = np.sum(Kp*(P-P0)**2)
    return E


## Different force models          
            
            
def force_vor_elastic(vor_regions,vor_point_region, vor_ridges,Kp,r0, cC, cV):
    # Calculates the force in all the vertices and points of the network:
    #vor is a dictionary with different values apparently or it can be a list of lists or something
    
    LC = len(cC)
    LV = len(cV)
    F_C = [] # Vector of resulting forces acting on the centers
    F_V = [] # Vector of resulting forces acting on the vertices
    
    dist_v_list = []
    
    
    for c in range(LC):
        Neigh_c = find_center_region(vor_regions,vor_point_region, c)
        
        f_c = 0.0*np.array([1.,1.])
        for n_c in Neigh_c:
            v_c = cV[n_c]
            if any(v_c < 5) and any(v_c > -5):
                r_vec_center = cC[c] - v_c
                abs_rvc = np.abs(r_vec_center)
                n_rvc = r_vec_center/abs_rvc
                f_c += -Kp*(abs_rvc-r0)*n_rvc
        F_C.append(f_c)
    
    for v in range(LV):
        f_v = 0.0*np.array([1.,1.])
        if any(cV[v] < 5) and any(cV[v] > -5): 
            NeighC_V, NeighV_V = find_vertex_neighbour(vor_regions,vor_point_region,vor_ridges,v)
            
            for nc_v in NeighC_V:
                c_v = cC[nc_v]
                r_center_vertex = cV[v] - c_v
                abs_rcv = np.abs(r_center_vertex)
                n_rcv = r_center_vertex/abs_rcv
            
                f_v += -Kp*(abs_rcv - r0)*n_rcv
        
            for nv_v in NeighV_V:
                v_v = cV[nv_v]
                if any(v_v < 5) and any(v_v > -5):
                    edgeV = cV[v] - v_v
                    abs_vv = np.abs(edgeV)
                    
                    dist_v_list.append(abs_vv)
                    n_vv = edgeV/abs_vv
            
                    f_v += -Kp*(abs_vv - r0)*n_vv
            
        F_V.append(f_v)
        
    return np.array(F_C), np.array(F_V), np.mean(dist_v_list)      


def phi_n(l_n):
    return l_n*np.tan(np.pi/l_n)


## T1 Transitions

def T1transition(vor_vertices, vor_ridges, vor_regions, vor_point_region,thresh_len):
    
    ''' This function runs through all the interior vertices on the network and does topological rearrangements (T1 transitions) if the length of the edges are below a certain threshold. The overall effect on the network should be similar to an viscoelastic relaxation '''
    
    #make a list of vertices to run through in the next cycle
    
    vertex_list = []
    
    for k in range(len(vor_vertices)):
        vertex_list.append(k)
        
    
    for i in vertex_list:
        
        if np.isnan(np.sum(vor_vertices[i])):
            continue
        else:
        
        
            #Find all neighbouring vertices of vertex i
        
            vertices_neigh = find_vertex_neighbour_vertices(vor_ridges,i)
        
            # list of neighbouring vertices lengths
            list_len = []
            list_neigh_v_not_excluded = []
            for v in vertices_neigh:
                if v in vertex_list:
                    list_neigh_v_not_excluded.append(int(v))
                    deltax = vor_vertices[i]-vor_vertices[int(v)]
                    list_len.append(np.sqrt(np.sum(deltax**2)))
        
            if len(list_neigh_v_not_excluded) > 2:
            
                # Find closest neighbouring vertices
        
                loc_v_min = np.argmin(list_len)
                lv_min = list_len[loc_v_min]
                v_min = int(list_neigh_v_not_excluded[loc_v_min])
        
                # Find neighbours of closest vertex to vertex i
        
                vertices_neigh_min = find_vertex_neighbour_vertices(vor_ridges,v_min)
            
                if len(vertices_neigh_min) > 2:
        
                    #Do rotation of the vertex edge along 90 degrees
                    th = np.pi/2
                    R_p = np.array([[0,-1],[1,0]])
                    R_m = np.array([[0,1],[-1,0]])
        
        
                    new_diff_v_vm_coord_1 = 0
                    new_diff_v_vm_coord_2 = 0
                    
        
                    if lv_min < thresh_len:
                        diff_v_vm = vor_vertices[i]-vor_vertices[v_min]
                        new_diff_v_vm_coord_1 = vor_vertices[v_min] + np.sqrt(3)/2*(R_p.dot(diff_v_vm)+1/2*diff_v_vm)
                        new_diff_v_vm_coord_2 = vor_vertices[v_min] + np.sqrt(3)/2*(R_m.dot(diff_v_vm)+1/2*diff_v_vm)
                        
                        with open('vertextransitions.txt','a') as text:
                    
                            text.write('Old length: ')
                    
                            text.write(str(lv_min)+' and the ')
                    
                            text.write('threshold: ')
                    
                            text.write(str(thresh_len)+'\n')
            
        
                        # Assign each vertex to each of the ends of the rotated vector
        
                        p_1 = np.random.rand()
        
                        if p_1 > 0.5:
                            vor_vertices[i] = new_diff_v_vm_coord_1
                            vor_vertices[v_min] = new_diff_v_vm_coord_2
                        elif p_1 <= 0.5:
                            vor_vertices[i] = new_diff_v_vm_coord_2
                            vor_vertices[v_min] = new_diff_v_vm_coord_1
            
        
                        #list of neighbouring vertice lengths in the closest neighbouring vertex for rotated vector only assign vertices in vertex list
        
                        list_len_vm = []
                        new_neighvm = []
                        
                        for v in vertices_neigh_min:
                            if v != i:
                                #For this one we don't need to remove previous vertices
                                deltax = vor_vertices[i]-vor_vertices[int(v)]
                                list_len_vm.append(np.sqrt(np.sum(deltax**2)))
                                new_neighvm.append(int(v))
        
                        list_len_v = []
                        new_neighv = []
                        for v in list_neigh_v_not_excluded:
                            if v != v_min:
                                deltax = vor_vertices[v_min]-vor_vertices[v]
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
            
                        loc_ridges = np.where(vor_ridges == i)[0]
            
                        loc_neigh_not_vm = np.where(np.array(vertices_neigh)!=v_min)[0]

                        skip_parameter = int(0)
                        for j in range(len(loc_ridges)):
                
                            if v_min in vor_ridges[loc_ridges[j]]:
                                skip_parameter += int(1)
                                continue
                            else:
                                js = int(j-skip_parameter)
                                vor_ridges[loc_ridges[j]]= [vertices_neigh[loc_neigh_not_vm[js]],i]                    
                        loc_ridges = np.where(vor_ridges == v_min)[0]
            
                        loc_neigh_not_i = np.where(np.array(vertices_neigh_min)!= i)[0]
            
                        skip_parameter = int(0)
                        
                        for j in range(len(loc_ridges)):
                            if i in vor_ridges[loc_ridges[j]]:
                                skip_parameter+=int(1)
                                continue
                            else:
                                
                                js = int(j-skip_parameter)
                                vor_ridges[loc_ridges[j]]= [vertices_neigh_min[loc_neigh_not_i[js]],v_min]
                                
                        
                        #Oh god, I forgot about the connections between the vertices and the centers which also change
                        
                        #For vertex i
                        regions_neigh_v = []
                        for v in vertices_neigh:
                            region, center = find_vertex_neigbour_centers(vor_regions,vor_point_region,v)
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
                            if i in vor_regions[r]:
                                continue
                            else:
                                #print(i)
                                #print(v_min)
                                #print(vor_regions[r])
                                
                                if v_min in vor_regions[r]:
                                    vor_regions[r].append(i)
                                    vor_regions[r].remove(v_min)
                                
                                #print(vor_regions[r])
                                
                        #For neighbouring vertex v_min
                        regions_neigh_vmin = []
                        for v in vertices_neigh_min:
                            region, center = find_vertex_neigbour_centers(vor_regions,vor_point_region,v)
                            regions_neigh_vmin.append(region)
                            
                        region_vmin = []
                        
                        for k in range(len(regions_neigh_vmin)):
                            region_vmin.append(list(set(regions_neigh_vmin[k]).intersection(regions_neigh_vmin[(k+1)%len(regions_neigh_v)])))
                            
                        region_vmin_true = []
                        
                        for r in region_vmin:
                            region_vmin_true = list(set(region_vmin_true).union(r))
                            
                        for r in region_vmin_true:
                            if v_min in vor_regions[r]:
                                continue
                            else:
                                if i in vor_regions[r]:
                                    vor_regions[r].remove(i)
                                    vor_regions[r].append(v_min)
                                
                        
                        
                            
                        
                            
                                
                    vertex_list.remove(i)
            elif len(list_neigh_v_not_excluded) <= 2:
                continue

    return vor_ridges, vor_vertices, vor_regions