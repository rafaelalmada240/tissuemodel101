import numpy as np
import selfpropelledparticlevoronoi as sppv
import pathlib
import h5py as h5
import matplotlib.pyplot as plt
import findconnections as fc

main_path = pathlib.Path().absolute()
datafileloadname = input('Name of savefile to load: ')
datafileloadname = datafileloadname + '.h5'
data_set = h5.File(str(main_path) + '/' + datafileloadname,'r')
coords = np.array(data_set['centers'])
vorVertices = np.array(data_set['vertices'])
vorRidges = np.array(data_set['Edge connections'])
Fvertex = np.array(data_set['forces_vertices'])

Fnorm = np.sqrt(np.sum(Fvertex**2,1))

#Find vertices in boundary
Bound = fc.find_boundary_vertices(len(vorVertices),vorRidges)

x_left = -2
x_right = 2

loc_bound_left = np.where(vorVertices[Bound,0,0]<x_left)[0]
bound_left = np.array(Bound)[loc_bound_left[np.argsort(vorVertices[loc_bound_left,1,0])]]
loc_bound_right = np.where(vorVertices[Bound,0,0]>x_right)[0]
bound_right = np.array(Bound)[loc_bound_right[np.argsort(vorVertices[loc_bound_right,1,0])]]

#plt.plot(vorVertices[bound_left,1,0])

Fnorm_maxl = Fnorm[bound_left,:1000].T/np.max(Fnorm[bound_left,:1000],1)
Fnorm_maxr = Fnorm[bound_right,:1000].T/np.max(Fnorm[bound_right,:1000],1)

plt.figure(figsize=(12,4))

plt.subplot(121)
plt.plot(Fnorm_maxl[:,len(bound_left)-1])
plt.plot(Fnorm_maxl[:,len(bound_left)-2])
plt.plot(Fnorm_maxl[:,len(bound_left)-3])
plt.plot(Fnorm_maxl[:,len(bound_left)-4])
plt.plot(Fnorm_maxl[:,len(bound_left)-5])
plt.plot(Fnorm_maxl[:,len(bound_left)-6])
plt.title('Force normalized to maximum')
plt.xlabel('Time')
#plt.plot(Fnorm[bound_left[len(bound_left)-3],:].T/np.max(Fnorm[bound_left[len(bound_left)-3],:1000]))
plt.legend(vorVertices[bound_left[::-1],1,0])
plt.xscale("log")
#plt.yscale("log")

plt.subplot(122)
plt.plot(Fnorm_maxr[:,len(bound_right)-1])
plt.plot(Fnorm_maxr[:,len(bound_right)-2])
plt.plot(Fnorm_maxr[:,len(bound_right)-3])
plt.plot(Fnorm_maxr[:,len(bound_right)-4])
plt.plot(Fnorm_maxr[:,len(bound_right)-5])
plt.plot(Fnorm_maxr[:,len(bound_right)-6])
plt.title('Force normalized to maximum')
plt.xlabel('Time')
#plt.plot(Fnorm[bound_left[len(bound_left)-3],:].T/np.max(Fnorm[bound_left[len(bound_left)-3],:1000]))
plt.legend(vorVertices[bound_right[::-1],1,0])
plt.xscale("log")
plt.savefig('NormForceShear1.png',dpi=150)


#displ = np.sqrt(np.sum(np.diff(vorVertices[bound_left],2)**2,1))

#print(vorVertices[bound_left,:,0][:,:,None].shape)
displ = np.zeros(np.sqrt(np.sum(vorVertices[bound_left],1)).shape)

for i in range(vorVertices.shape[2]):
    displ[:,i] = np.sqrt(np.sum((vorVertices[bound_left,:,0]-vorVertices[bound_left,:,i])**2,1))
    
dispr = np.zeros(np.sqrt(np.sum(vorVertices[bound_right],1)).shape)

for i in range(vorVertices.shape[2]):
    dispr[:,i] = np.sqrt(np.sum((vorVertices[bound_right,:,0]-vorVertices[bound_right,:,i])**2,1))



plt.figure(figsize=(12,4))

plt.subplot(121)
plt.plot(displ[len(bound_left)-1,:1000].T,Fnorm_maxl[:,len(bound_left)-1],'.-')
plt.plot(displ[len(bound_left)-2,:1000].T,Fnorm_maxl[:,len(bound_left)-2],'.-')
plt.plot(displ[len(bound_left)-3,:1000].T,Fnorm_maxl[:,len(bound_left)-3],'.-')
plt.plot(displ[len(bound_left)-4,:1000].T,Fnorm_maxl[:,len(bound_left)-4],'.-')
plt.plot(displ[len(bound_left)-5,:1000].T,Fnorm_maxl[:,len(bound_left)-5],'.-')
plt.plot(displ[len(bound_left)-6,:1000].T,Fnorm_maxl[:,len(bound_left)-6],'.-')
plt.legend(vorVertices[bound_left[::-1],1,0])
plt.title('Force normalized to maximum')
plt.xlabel('displacement')
#plt.plot(displ[len(bound_left)-1,:1000].T,Fnorm[bound_left[len(bound_left)-1],:1000].T/np.max(Fnorm[bound_left[len(bound_left)-1],:1000]),'.-')
plt.xscale("log")
#plt.yscale("log")

plt.subplot(122)
plt.plot(dispr[len(bound_right)-1,:1000].T,Fnorm_maxr[:,len(bound_right)-1],'.-')
plt.plot(dispr[len(bound_right)-2,:1000].T,Fnorm_maxr[:,len(bound_right)-2],'.-')
plt.plot(dispr[len(bound_right)-3,:1000].T,Fnorm_maxr[:,len(bound_right)-3],'.-')
plt.plot(dispr[len(bound_right)-4,:1000].T,Fnorm_maxr[:,len(bound_right)-4],'.-')
plt.plot(dispr[len(bound_right)-5,:1000].T,Fnorm_maxr[:,len(bound_right)-5],'.-')
plt.plot(dispr[len(bound_right)-6,:1000].T,Fnorm_maxr[:,len(bound_right)-6],'.-')
plt.legend(vorVertices[bound_right[::-1],1,0])
plt.title('Force normalized to maximum')
plt.xlabel('displacement')
#plt.plot(displ[len(bound_left)-1,:1000].T,Fnorm[bound_left[len(bound_left)-1],:1000].T/np.max(Fnorm[bound_left[len(bound_left)-1],:1000]),'.-')
plt.xscale("log")
#plt.yscale("log")
plt.savefig('NormForceShearDisp1.png',dpi=150)


loc_Fmaxl = np.argmax(Fnorm_maxl,0)

Fmaxl = []

displ_loc = []

ylocl = []

for i in range(len(loc_Fmaxl)):
    Fmaxl.append(Fnorm[bound_left[i],loc_Fmaxl[i]])
    displ_loc.append(displ[i,loc_Fmaxl[i]])
    ylocl.append(vorVertices[bound_left[i],1,0])
    
loc_Fmaxr = np.argmax(Fnorm_maxr,0)

Fmaxr = []

dispr_loc = []

ylocr = []

for i in range(len(loc_Fmaxr)):
    Fmaxr.append(Fnorm[bound_right[i],loc_Fmaxr[i]])
    dispr_loc.append(dispr[i,loc_Fmaxr[i]])
    ylocr.append(vorVertices[bound_right[i],1,0])


plt.figure(figsize=(20,4))

plt.subplot(141)
plt.plot(np.array(ylocl)[7:],np.array(Fmaxl)[7:],'o--')
plt.plot(np.array(ylocr)[7:],np.array(Fmaxr)[7:],'o--')
plt.title('Maximum Force vs height')
plt.yscale("log")

plt.subplot(142)
plt.plot(np.array(ylocl)[7:],np.array(displ_loc)[7:],'o--')
plt.plot(np.array(ylocr)[7:],np.array(dispr_loc)[7:],'o--')
plt.title('displacement at maximum force vs height')
#plt.yscale("log")

plt.subplot(143)
plt.plot(loc_Fmaxl[7:],np.array(Fmaxl)[7:],'o--')
plt.plot(loc_Fmaxr[7:],np.array(Fmaxr)[7:],'o--')
plt.title('Maximum Force vs time')
plt.yscale("log")

plt.subplot(144)
plt.plot(loc_Fmaxl[7:],np.array(displ_loc)[7:],'o--')
plt.plot(loc_Fmaxr[7:],np.array(dispr_loc)[7:],'o--')
plt.title('displacement at maximum force vs time')
#plt.plot(displ[len(bound_left)-3,:].T)
#plt.plot(displ[len(bound_left)-1,:].T)
#plt.xscale("log")

plt.savefig('MaxForceDisp.png',dpi=150)

loc_displmax = np.argmax(displ[:,:1000],1)

displ_max = []

for i in range(len(loc_displmax)):
    
    displ_max.append(displ[i,loc_displmax[i]])
    
loc_disprmax = np.argmax(dispr[:,:1000],1)

dispr_max = []

for i in range(len(loc_disprmax)):
    
    dispr_max.append(dispr[i,loc_disprmax[i]])



plt.figure()
plt.plot(np.array(ylocl)[7:],np.array(displ_max)[7:],'o--')
plt.plot(np.array(ylocr)[7:],np.array(dispr_max)[7:],'o--')
plt.title('Maximum displacement vs height')
plt.savefig('MaxDispHeight.png',dpi=150)

for i in range(len(ylocr)):
    with open('MaxDispHeight.txt','a') as text:
        text.write(str(ylocr[i]))
        text.write(' ')
        text.write(str(dispr_max[i]))
        text.write(' ')
        text.write(str(ylocl[i]))
        text.write(' ')
        text.write(str(displ_max[i])+ '\n')
#plt.yscale("log")



