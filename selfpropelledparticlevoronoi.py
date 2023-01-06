import numpy as np

def find_center_region(vor,center):
    #Gives all the neighbours of a center in a voronoi tesselation (removes all -1 vertices)
    R = vor.regions[vor.point_region[center]]
    for e in R:
        if e==-1:
            R.remove(e)
    return R



def find_vertex_neighbour(vor,vertex):
    
    #Gives all the neighbours of a vertex in a voronoi tesselation (removes all -1 vertices)
    
    list_regions = []
    list_vertex_neigh = []
    list_centers = []
    i = 0
    for region in vor.regions:
        if vertex in region:
            list_regions.append(i)
        i+=1       
    for region in list_regions:
        loc_points = np.where(vor.point_region==region)
        list_centers.append(loc_points[0][0])
        
    for ridge in vor.ridge_vertices:
        if vertex in ridge:            
            for elem in ridge:
                if (elem != vertex) and (elem != -1):
                    list_vertex_neigh.append(elem)
    
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
            
            
def force_vor_elastic(vor,Kp,r0, cC, cV):
    # Calculates the force in all the vertices and points of the network:
    
    C = vor.points
    V = vor.vertices
    F_C = [] # Vector of resulting forces acting on the centers
    F_V = [] # Vector of resulting forces acting on the vertices
    
    
    for c in range(len(C)):
        Neigh_c = find_center_region(vor, c)
        
        f_c = 0.0*np.array([1.,1.])
        for n_c in Neigh_c:
            v_c = cV[n_c]
            if any(v_c < 5) and any(v_c > -5):
                r_vec_center = cC[c] - v_c
                abs_rvc = np.abs(r_vec_center)
                n_rvc = r_vec_center/abs_rvc
                f_c += -Kp*(abs_rvc-r0)*n_rvc
        F_C.append(f_c)
    
    for v in range(len(V)):
        f_v = 0.0*np.array([1.,1.])
        if any(cV[v] < 5) and any(cV[v] > -5): 
            NeighC_V, NeighV_V = find_vertex_neighbour(vor,v)
            
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
                    n_vv = edgeV/abs_vv
            
                    f_v += -Kp*(abs_vv - r0)*n_vv
            
        F_V.append(f_v)
        
    return np.array(F_C), np.array(F_V)      


def phi_n(l_n):
    return l_n*np.tan(np.pi/l_n)

