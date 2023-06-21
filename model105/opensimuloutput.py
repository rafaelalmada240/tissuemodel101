import numpy as np
import matplotlib.pyplot as plt

p_list = []
a_list = []
for lr in range(37):
    for lw in range(10):
        with open('simul/woundinfoA2.900469900595457G1L-'+str(lr/10)+'Lw-'+str(lw)+'Ncells225.txt','r') as text:
            for line in text:
                pass
            last_line = (line.replace("\n","")).split(' ')
            p_list.append([lr/10,lw,float(last_line[1])])
            a_list.append([lr/10,lw,float(last_line[2])])

plt.figure(figsize=(16,4))
a_array = np.array(a_list)
for i in range(185):
    if a_array[i,2]>2:
        if i == 0:
            plt.plot(a_array[i,0],a_array[i,1],'go',alpha=0.6,label = "wound stays open")
        else:
            plt.plot(a_array[i,0],a_array[i,1],'go',alpha=0.6)
    if a_array[i,2]<2:
        if i == 184:
            plt.plot(a_array[i,0],a_array[i,1],'rd',alpha=0.6,label = "wound closes")
        else:
            plt.plot(a_array[i,0],a_array[i,1],'rd',alpha=0.6)
plt.xlabel("$\\beta$",fontsize = 18)
plt.ylabel("$\lambda_W$",fontsize = 18)
plt.title("Wound area",fontsize=24)
plt.legend(fontsize=16)


a_array2 = np.zeros((37,5))
for i in range(a_array.shape[0]):
    a_array2[i//5,i%5] = a_array[i,2]

plt.plot(a_array2[:,0],'.-',label="$\lambda_W = 0$")
plt.plot(a_array2[:,1],'o-',label="$\lambda_W = 1$")
plt.plot(a_array2[:,2],'s-',label="$\lambda_W = 2$")
plt.plot(a_array2[:,3],'d-',label="$\lambda_W = 3$")
plt.plot(a_array2[:,4],'-',label="$\lambda_W = 4$")
plt.legend(fontsize=16)
plt.xlabel("$10\\beta$",fontsize = 18)
