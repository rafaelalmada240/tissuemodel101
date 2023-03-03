import numpy as np
import selfpropelledparticlevoronoi as sppv


## Topological properties

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

def graph_energy(Amatrix):
    eig_spectrum = np.linalg.eig(Amatrix)[0]
    return np.sum(np.abs(eig_spectrum))
    

## T1 Transitions

def line(a,b,x):
    return a*x+b

def line_coeffs(p1,p2):
    a = (p2[1]-p1[1])/(p2[0]-p1[0])
    b = p1[1]-a*p1[0]
    return a,b

def T1_rotations(vertices, v_orig,v_min):
    '''Rotates the coordinates of an edge undergoing a T1 transition
    Input:
    vertices - set of coordinate values for each vertex
    (v_orig,v_min) - edge in the iteration that is undergoing a T1 transition
    
    Output:
    (r_vorig, r_vmin) - new coordinate values of the rotated edge
    '''
    th = np.pi/2
    
    old_vorig = vertices[v_orig]
    old_vmin = vertices[v_min]
    
    a1, b1 = line_coeffs(old_vorig,old_vmin)
    R_p = np.array([[0,-1],[1,0]])
    R_m = np.array([[0,1],[-1,0]])
                    
    diff_v_vm = old_vorig-old_vmin
    new_diff_v_vm_coord_1 = old_vmin + np.sqrt(3)/2*(R_p.dot(diff_v_vm)+1/2*diff_v_vm)
    new_diff_v_vm_coord_2 = old_vmin + np.sqrt(3)/2*(R_m.dot(diff_v_vm)+1/2*diff_v_vm)
            
        
    # Assign each vertex to each of the ends of the rotated vector
        
    p_1 = np.random.rand()
     
    r_vorig = old_vorig
    r_vmin = old_vmin 
        
    if p_1 > 0.5:
        r_vorig = new_diff_v_vm_coord_1
        r_vmin = new_diff_v_vm_coord_2
    else:
        r_vorig = new_diff_v_vm_coord_2
        r_vmin = new_diff_v_vm_coord_1
        
    a2, b2 = line_coeffs(r_vorig,r_vmin)
        
    return r_vorig, r_vmin, (a1,b1),(a2,b2), old_vorig, old_vmin

def T1_orientation_matrix(v_aold,v_bold, v_anew,v_bnew):
    
    #In theory the midpoints should be equal so the average is unnecessary, but it's better to be safe than sorry
    mid_point = 0.5*(v_bold+0.5*(v_aold - v_bold))+0.5*(v_bnew+0.5*(v_anew - v_bnew))
    
    raold = v_aold - mid_point
    ranew = v_anew - mid_point
    
    Ainv = np.array([raold,ranew]).T
    
    A = np.linalg.inv(Ainv)
    
    return A, mid_point

def T1_len_neigh_after_rot(vertices, v_orig,neigh_min,v_min,neigh_orig, coefs1,coefs2,rvorig0, rvmin0):
    '''
    During a T1 transition finds the vertices for the new neighbourhood for the edge undergoing a T1 transition, using projection (in theory should be robust to several configurations)
    Input:
    vertices - set of coordinate values for each vertex
    (v_orig,v_min) - edge in the iteration that is undergoing a T1 transition
    neigh_min - neighbourhood of v_min
    neigh_orig - neighbourhood of v_orig
    
    Output:
    new_neighv,new_neighm - Neighbourhoods of vorig and vmin that are not vorig or vmin, respectively
    lvv_min, lvm_min - the index of the vertices in new_neighv and new_neighm that will be swapped between the two neighbourhoods, respectively
    '''
    
    TransfMat, midpoint = T1_orientation_matrix(rvorig0,rvmin0,vertices[v_orig],vertices[v_min])
    
    va0 = TransfMat.dot(rvorig0-midpoint)
    vb0 = TransfMat.dot(rvmin0-midpoint)
    
    va1 = TransfMat.dot(vertices[v_orig]-midpoint)
    vb1 = TransfMat.dot(vertices[v_min]-midpoint)
    
    #Surprisingly longer than expected, will also have to define up and down
    
    f1_neighm = []
    
    list_m_len = [] #This should go from -1 to 1
    neighm_ind = []
                
    for v in neigh_min:
        if v != v_orig:
            r_neighm = TransfMat.dot(vertices[int(v)]-midpoint)
            lr_check = int(va1.dot(r_neighm-vb0)>0)
            ud_check = int(va0.dot(r_neighm-vb0)>0)
            list_m_len.append(va0.dot(r_neighm-vb0)/((r_neighm-vb0).dot(r_neighm-vb0))**0.5)
            f1_neighm.append(lr_check)
            
            neighm_ind.append(int(v))
     
   
    if np.diff(np.array(f1_neighm))==0:
        #A bunch of conditions
        lm_min = np.argmin(np.array(list_m_len))
        f1_neighm[lm_min]=1-f1_neighm[lm_min] 
       
    f1_neighv = []
    f2_neighv = []
    list_v_len = [] #This should go from -1 to 1
    neighv_ind = []
                
    for v in neigh_orig:
        if v != v_min:
            r_neighv = TransfMat.dot(vertices[int(v)]-midpoint)
            lr_check = int(va1.dot(r_neighv-va0)>0)
            ud_check = int(va0.dot(r_neighv-va0)>0)
            list_v_len.append(va0.dot(r_neighv-va0)/((r_neighv-va0).dot(r_neighv-va0))**0.5)
            f1_neighv.append(lr_check)
            f2_neighv.append(ud_check)
            
            neighv_ind.append(int(v))
            
     #In case of weird freak configurations (both in same side LR)
     
    if np.diff(np.array(f1_neighv))==0:
        lv_min = np.argmax(np.array(list_v_len))
        f1_neighv[lv_min]=1-f1_neighv[lv_min]  
 
    lvm_min = np.where(np.array(f1_neighm)==0)[0][0]
    lvv_min = np.where(np.array(f1_neighv)==1)[0][0]

    return neighv_ind,neighm_ind, lvv_min, lvm_min

def T1_change_edges(ridges,vertices_neigh, vertices_neigh_min,i,v_min):
    
    '''
    After the change in neighbours for each vertex in the edge that was rearranged, remove the old edges from the edge set and include the new edges
    '''
    
    
    loc_ridges = np.where(ridges == i)[0]
    
    #print(loc_ridges)
            
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
    #print((len(loc_ridges),len(loc_neigh_not_i),len(vertices_neigh_min)))
            
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
        
            vertices_neigh = sppv.find_vertex_neighbour_vertices(ridges,i)
        
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
        
                vertices_neigh_min = sppv.find_vertex_neighbour_vertices(ridges,v_min)
                
                #Only do this for vertices with 3 neighbours, and also avoid triangles
            
                if (len(vertices_neigh_min) > 2) and (len(list(set(list_neigh_v_not_excluded).intersection(vertices_neigh_min)))<1):
                    
                    #Do rotation of the vertex edge along 90 degrees (The geometric part is okay, not likely to change much)
                    if lv_min < thresh_len: 
                        vertices[i], vertices[v_min], coefs1, coefs2, v_oldi, vold_min = T1_rotations(vertices,i,v_min)
                        
                        #list of neighbouring vertice lengths in the closest neighbouring vertex for rotated vector only assign vertices in vertex list        
                        new_neighv, new_neighvm, lvv_min, lvm_min = T1_len_neigh_after_rot(vertices,i,vertices_neigh_min,v_min,list_neigh_v_not_excluded, coefs1,coefs2, v_oldi, vold_min)   
                               
                        vertices_neigh_min = list(set(vertices_neigh_min).difference({new_neighvm[lvm_min]}).union({new_neighv[lvv_min]}))
                        vertices_neigh = list(set(list_neigh_v_not_excluded).difference({new_neighv[lvv_min]}).union({new_neighvm[lvm_min]}))
                               
                        #Wait, you actually need to change the ridges as well, hopefully this works                        
                        ridges = T1_change_edges(ridges,vertices_neigh,vertices_neigh_min,i,v_min)
                                   
                        #Oh god, I forgot about the connections between the vertices and the centers which also change                       
                        #For vertex i
                        regions_neigh_v = []
                        for v in vertices_neigh:
                            ff = sppv.find_vertex_neighbour_centers(regions,point_region,v)
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
                            ff = sppv.find_vertex_neighbour_centers(regions,point_region,v)
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
