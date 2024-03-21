import numpy as np
import matplotlib.pyplot as plt
import movies_from_plots as mfp
from scipy.spatial import Voronoi, voronoi_plot_2d
import selfpropelledparticlevoronoi as sppv
import findconnections as fc
#import pathlib
#import h5py as h5
import psutil as ps
from IPython.core.display import HTML
import time

#from https://www.kaggle.com/getting-started/210022 with some changes
def restart_kernel_and_run_all_cells():
    HTML(
        '''
            <script>
                code_show = false;
                function restart_run_all(){
                    IPython.notebook.kernel.restart();
                }
                function code_toggle() {
                    if (code_show) {
                        $('div.input').hide(200);
                    } else {
                        $('div.input').show(200);
                    }
                    code_show = !code_show
                }
                code_toggle() 
                restart_run_all()
            </script>

        '''
    )
    

# Plot frames from different time steps   


def open_file(foldername0,N,issubfolder):
    
    coords_list = []
    point_regions = []
    regions_list = []
    vertex_list = []
    boundaries_list = []
    edges_list = []
    for i in range(N):
        #print(i)
        coords = []
        vorPointRegion1 = []
        vorRegions = []
        
        if issubfolder == 0:
            foldername = foldername0+'movie_output/'
        else:
            foldername = foldername0

        with open(foldername+'centers'+str(i)+'.txt','r') as text:
                    for line in text:
                        last_line = (line.replace("\n","")).split(';')
                        coords.append([float(last_line[1]),float(last_line[2])])
                        vorPointRegion1.append(int(last_line[3]))
                        l4 = last_line[4].replace("[","").replace("]","").split(',')
                        lint4 = []
                        if all('' == s or s.isspace() for s in l4):
                            vorRegions.append(lint4)
                            continue
                        else:
                            for r in l4:
                                lint4.append(int(r))
                        vorRegions.append(lint4)
                        
        coords = np.array(coords)

        vertices = []
        with open(foldername+'vertices'+str(i)+'.txt','r') as text:
                    for line in text:
                        last_line = (line.replace("\n","")).split(';')
                        vertices.append([float(last_line[1]),float(last_line[2])])
        vertices = np.array(vertices)

        #print(vorPointRegion)

        vorRidges = []
        with open(foldername+'edges'+str(i)+'.txt','r') as text:
                    for line in text:
                        last_line = (line.replace("\n","")).split(';')
                        l4 = last_line[1].replace("[","").replace("]","").split(',')
                        lint4 = []
                        if all('' == s or s.isspace() for s in l4):
                            vorRidges.append(lint4)
                            continue
                        else:
                            for r in l4:
                                lint4.append(int(r))
                        vorRidges.append(np.array(lint4))
                        

        wloc = 0

        with open(foldername+'woundloc'+str(i)+'.txt','r') as text:
            for line in text:
                wloc = int(line.replace("\n",""))

        vorPointRegion= []
        for k in range(len(coords)):
            vorPointRegion.append(k)

        Boundaries = []                
        with open(foldername+'boundaries'+str(i)+'.txt','r') as text:
            for line in text:
                last_line = line.replace("\n","")
                l4 = line.replace("[","").replace("]","").split(',')
                lint4 = []
                if all('' == s or s.isspace() for s in l4):
                    Boundaries.append(lint4)
                    continue
                else:
                    for r in l4:
                        lint4.append(int(r))
                Boundaries.append(lint4)
                
        coords_list.append(coords)
        point_regions.append(vorPointRegion)
        regions_list.append(vorRegions)
        vertex_list.append(vertices)
        boundaries_list.append(Boundaries)
        edges_list.append(vorRidges)
    dataset = {"centers": coords_list, "vertices": vertex_list, "Edge connections":edges_list
               , "WoundLoc":wloc,"boundaries":boundaries_list,"regions":regions_list,"point regions": point_regions}
    return dataset


#Add regions and wloc to the mix:


def saveframe(x_array, v_array, ridge_array_1, boundaries, N_init,N_fin,f_step,file_path,regions, wloc):
    t_initial = time.time() #time where the printing of frames starts
    fig = plt.figure(figsize=(8,8))
    #ordered_boundary1 = sppv.rearrange(len(boundaries[0][1]),sppv.adj_mat(boundaries[0][1],ridge_array))
    #ordered_boundary1.append(ordered_boundary[0])
    #plt.plot(v_array[0,list(np.array(boundaries[0][1])[ordered_boundary]),0],v_array[0,list(np.array(boundaries[0][1])[ordered_boundary]),1],'--', color='k', lw = 3)  
    for k in range(N_init,N_fin,f_step):
        ridge_array = ridge_array_1[k,:,:].astype(int)
        
        
        V = Voronoi(x_array[k,:,:])
        # voronoi_plot_2d(V,show_vertices=False,show_points=False,line_alpha=0.1)
        ordered_boundary = sppv.rearrange(len(boundaries[k][1]),sppv.adj_mat(boundaries[k][1],ridge_array))
        ordered_boundary.append(ordered_boundary[0])
        for i in range(len(ridge_array)):
           
            if any(np.array(ridge_array[i])==-1)==0:
                plt.plot(v_array[k,ridge_array[i],0],v_array[k,ridge_array[i],1],'-',color='lawngreen',alpha=1,lw=3)
        # plt.plot(v_array[k,:,0],v_array[k,:,1],'r.',alpha=0.5)
        
        plt.plot(v_array[k,list(np.array(boundaries[k][1])[ordered_boundary]),0],v_array[k,list(np.array(boundaries[k][1])[ordered_boundary]),1],'-', color='fuchsia', lw = 5)  
        #
        polygon = [[-5,-5],[-5,5],[5,5],[5,-5]]
        plt.fill(*zip(*polygon),'midnightblue',alpha=0.75)
        # for region in regions[k]:
        # 	if not -1 in region:	
        # 		polygon = [v_array[k,j] for j in region]
        # 		if len(polygon) <= 3:
        # 			plt.fill(*zip(*polygon),'darkviolet',alpha=0.5)
        # 		if len(polygon) == 4:
        # 			plt.fill(*zip(*polygon),'midnightblue',alpha=0.5)
        # 			#plt.fill(*zip(*polygon),'midnightblue',alpha=0.5)
        # 		if len(polygon) == 5:
        # 			plt.fill(*zip(*polygon),'deepskyblue',alpha=0.5)
        # 			#plt.fill(*zip(*polygon),'deepskyblue',alpha=0.5)
        # 		if len(polygon) == 6:
        # 			plt.fill(*zip(*polygon),'forestgreen',alpha=0.5)
        # 			#plt.fill(*zip(*polygon),'forestgreen',alpha=0.5)
        # 		if len(polygon) >= 7:
        # 			plt.fill(*zip(*polygon),'darkorange',alpha=0.5)
        # 			#plt.fill(*zip(*polygon),'darkorange',alpha=0.5)
        
        wound_poly =  [v_array[k,i] for i in list(np.array(boundaries[k][1])[ordered_boundary]) if fc.norm(v_array[k,i]) < np.sqrt(2)*5]    
        plt.fill(*zip(*wound_poly),'r',alpha=0.75)
        plt.figtext(0.5,0.8,'$\mathregular{\lambda_W = 5.}$',fontsize = 50,color='white')
        plt.figtext(0.5,0.7,'$\mathregular{p_0 = 4.}$',fontsize = 50,color='white')
        
        #plt.plot(x_array[k,:,0],x_array[k,:,1],'b.')
        
        

        # plt.xlabel('x')
        # plt.ylabel('y')
        # plt.title('Step = '+str(k))
        
        plt.xlim(-4,4)
        plt.ylim(-4,4)
        # plt.grid('False')
        ax = plt.gca()
        ax.set_xticks([])
        ax.set_yticks([])   
        plt.savefig(file_path+str(int(k/f_step))+'.png',dpi=150,bbox_inches='tight')
        fig.clear()
    
        
        
        
        memory_val = ps.virtual_memory()
        memory_percent=memory_val.percent

        memory_thresh = 50.0
        if memory_percent>=memory_thresh:
            print('Memory overload')
            t_final = time.time() - t_initial
            with open('log.txt','w') as log:
                log.write('Memory overloaded: Simulation stops at: '+str(t_final)+'s \n')
                log.write('frame where memory crashed: '+str(k))
            restart_kernel_and_run_all_cells()
            #HTML("<script>Jupyter.notebook.kernel.restart()</script>")
            break
    plt.close(fig)
    #plt.show()
    return

def makemovie(x_array, v_array, ridge_array,boundaries,regions,wloc):
    N_initial = int(input('Starting frame: '))
    N_final = int(input('Stopping frame: '))
    fr_step = int(input('Frame step: '))
    file_path = 'f'
    
    saveframe(x_array, v_array, ridge_array, boundaries, N_initial,N_final,fr_step,file_path,regions, wloc)
    return

# Load dataset to use

#main_path = pathlib.Path().absolute()
datafileloadname = input('Number of points to open: ')
#datafileloadname = datafileloadname + '.h5'
# foldername1 = input('Where is the movieoutput?: ')
foldername2 = input('Which tissue?: ')
foldername3 = input('Which wound size?: ')
issubfolder = int(input('Is it in a sub-folder (y-1,n-0): '))
if issubfolder==1:
    whichl = input('p0? - ')
    whichlw = input('lw? - ')
    foldername1 = 'l'+whichl+'lw'+whichlw+'/'
    foldername = 'tissues/tissue'+foldername2+'/size'+foldername3+'/movie_output/'+foldername1
else:
    foldername = 'tissues/tissue'+foldername2+'/size'+foldername3+'/'

data_set = open_file(foldername,int(datafileloadname),issubfolder)

coords_evo = np.array(data_set['centers'])
print(coords_evo.shape)
coords_evo_vertex = np.array(data_set['vertices'])
print(coords_evo_vertex.shape)
ridge_vectors = np.array(data_set['Edge connections'])
print(ridge_vectors.shape)
print(ridge_vectors[0,1])
boundary = data_set['boundaries']
print(boundary[1])

#Additional code for polygon coloring
wloc = data_set['WoundLoc']
regions = data_set['regions']

# Make the movie
makemovie(coords_evo,coords_evo_vertex, ridge_vectors,boundary,regions, wloc)

prompt = int(input('Do you want to make the movie now [y(1),n(0)]: '))
if prompt == 1:
        
    N_init = 0
    N_end = int(input('Stopping frame: '))
    fr_step = 1
    file_path = 'f'
    filename = input('Video file name (with .mp4 included): ')
        
    img_array,size = mfp.loadframe(N_init,N_end,fr_step,file_path)
    mfp.savevideo(img_array,filename,size,file_path)
