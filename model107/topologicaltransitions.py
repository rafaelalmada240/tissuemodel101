import numpy as np
import selfpropelledparticlevoronoi as sppv
from scipy.spatial import Voronoi
import findconnections as fc

## Topological properties




    

## T1 Transitions

def T1_change_edges(ridges,vertices_neigh, vertices_neigh_min,i,v_min):
    
    '''
    After the change in neighbours for each vertex in the edge that was rearranged, remove the old edges from the edge set and include the new edges
    '''
    
    
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
    




#To Change (switch cells of new configuration before switching vertices)

def T1transition2(vertices, ridges, regions, point_region,thresh_len):
    
    ''' This function runs through all the interior vertices on the network and does topological rearrangements (T1 transitions) 
    Using a set operation approach (change regions before doing the edge swap)
    
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
        
            vertices_neigh = fc.find_vertex_neighbour_vertices(ridges,i)
        
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
        
                vertices_neigh_min = fc.find_vertex_neighbour_vertices(ridges,v_min)
                
                #Only do this for vertices with 3 neighbours, and also avoid triangles
            
                if (len(vertices_neigh_min) > 2) and (len(list(set(list_neigh_v_not_excluded).intersection(vertices_neigh_min)))<1):
                    
                    #For vertex i
                    regions_neigh_v, center = fc.find_vertex_neighbour_centers(regions,point_region,i)
                                
                    #For neighbouring vertex v_min
                    regions_neigh_vmin, center = fc.find_vertex_neighbour_centers(regions,point_region,v_min)

                    #Find neighbouring regions prior to transition
                        
                    region_common = list(set(regions_neigh_v).intersection(regions_neigh_vmin))
                        
                    region_exc_v = list(set(regions_neigh_v).difference(region_common))
                
                    region_exc_vmin = list(set(regions_neigh_vmin).difference(region_common))
                    
                    if ((len(region_exc_v)>0) and (len(region_exc_vmin)>0)) and (len(region_common)>1): 
                        if lv_min < thresh_len: 
                            #Oh god, I forgot about the connections between the vertices and the centers which also change      
        
                        
                            #Topological rearrangement is equivalent to change of neighbours.
                            new_region_common = list(set(region_exc_v).union(region_exc_vmin))

                            p_1 = np.random.rand()
                            i_v = 0
                            i_min = 1
                            if p_1 < 0.5:
                                i_v = 1
                                i_min = 0
                        
                            #New exclusive regions for i and v_min

                            
                            new_region_exc_v = region_common[i_v]   
                            new_region_exc_min = region_common[i_min]

                        
                            new_region_v = list(set(new_region_common).union([new_region_exc_v]))      
                            new_region_vmin = list(set(new_region_common).union([new_region_exc_min]))
                        
                            #Removing and adding elements in new regions list

                            regions[region_exc_vmin[0]].append(i)
                            regions[region_exc_v[0]].append(v_min)
                        
                            
                        
                            regions[new_region_exc_v].remove(v_min)
                            regions[new_region_exc_min].remove(i)
                            

                            
                            
                            edge_common_v = list(set(regions[new_region_common[0]]).intersection(regions[new_region_common[1]]))

                        
                            new_neighv = list(set(regions[new_region_exc_v]).intersection(regions[region_exc_vmin[0]]).difference(set(edge_common_v).intersection(regions[new_region_exc_v])))

                            new_neighvm = list(set(regions[region_exc_v[0]]).intersection(regions[new_region_exc_min]).difference(set(edge_common_v).intersection(regions[new_region_exc_min])))
                            

                            #list of neighbouring vertice lengths in the closest neighbouring vertex for rotated vector only assign vertices in vertex list  

                            vertices_neigh_min = list(set(vertices_neigh_min).difference(new_neighv).union(new_neighvm))
                            vertices_neigh = list(set(vertices_neigh).difference(new_neighvm).union(new_neighv))
                            


                            #Wait, you actually need to change the ridges as well, I forgot it                        
                            ridges1 = T1_change_edges(ridges,vertices_neigh,vertices_neigh_min,i,v_min) 
                        
                            transition_counter += 1
                        
                            vertex_list.remove(i)
                            list_ridges1 = list(ridges1)
                            ridgesl = [list(l) for l in list_ridges1]
                        
                    
            elif len(list_neigh_v_not_excluded) <= 2:
                ridgesl =[list(l) for l in list(ridges)] 
                continue

    return ridgesl, vertices, regions, transition_counter