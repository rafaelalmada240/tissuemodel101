import numpy as np
import matplotlib.pyplot as plt
import movies_from_plots as mfp
import pathlib
import h5py as h5
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


def saveframe(x_array, v_array, ridge_array, N_init,N_fin,f_step,file_path):
    t_initial = time.time() #time where the printing of frames starts
    for k in range(N_init,N_fin,f_step):
        plt.figure(figsize=(16,16))

        for i in range(len(ridge_array)):
           
            if any(np.array(ridge_array[i])==-1)==0:
                plt.plot(v_array[ridge_array[i],0,k],v_array[ridge_array[i],1,k],'g-',alpha=0.3)
        plt.plot(v_array[:,0,k],v_array[:,1,k],'r.',alpha=0.5)
        plt.plot(x_array[:,0,k],x_array[:,1,k],'b.')

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Step = '+str(k))
        
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.grid('both')
        plt.savefig(file_path+str(int(k/f_step))+'.png',dpi=100)
        plt.close()
        
        
        
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
    return

def makemovie(x_array, v_array, ridge_array):
    N_initial = int(input('Starting frame: '))
    N_final = int(input('Stopping frame: '))
    fr_step = int(input('Frame step: '))
    file_path = input('Image File Path: ')
    
    saveframe(x_array, v_array, ridge_array, N_initial,N_final,fr_step,file_path)
    return

# Load dataset to use

main_path = pathlib.Path().absolute()
datafileloadname = input('Name of savefile: ')
datafileloadname = datafileloadname + '.h5'
data_set = h5.File(str(main_path) + '/' + datafileloadname,'r')

coords_evo = np.array(data_set['centers'])
coords_evo_vertex = np.array(data_set['vertices'])
ridge_vectors = np.array(data_set['Edge connections'])

# Make the movie
makemovie(coords_evo,coords_evo_vertex, ridge_vectors)

prompt = int(input('Do you want to make the movie now [y(1),n(0)]: '))
if prompt == 1:
        
    N_init = int(input('Starting frame: '))
    N_end = int(input('Stopping frame: '))
    fr_step = int(input('Frame step: '))
    file_path = input('Image File Path: ')
    filename = input('Video file name (with .mp4 included): ')
        
    img_array,size = mfp.loadframe(N_init,N_end,fr_step,file_path)
    mfp.savevideo(img_array,filename,size,file_path)