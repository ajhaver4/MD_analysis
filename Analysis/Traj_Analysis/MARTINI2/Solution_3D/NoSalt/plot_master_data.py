import pandas
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
file=sys.argv[1]
data = pandas.read_csv(file,comment='#',sep='\t')


n_values = len(data['d1'])
timesteps_array = data['Timestep']/1e6




# print("lenght of frames: ",len(min_arr))
# sys.exit()
# timesteps_array = data['Frame No.']

"""
Uncomment for selecting a time range of values
"""
#Plotting CVs d1 and d2
# mask = (data['Timestep'] <129000000) & (data['Min dist']<0.6)
# n_values = len(data['d1'])
# timesteps_array = data['Timestep'][mask]

# fig_dist,ax_dist = plt.subplots(2,1,sharex=True)
# ax_dist[0].plot(timesteps_array,data['d1'][mask],'r',linewidth=0.6,label='d1')
# ax_dist[0].plot(timesteps_array,data['d2'][mask],'b',linewidth=0.6,label='d2')
# ax_dist[0].plot(timesteps_array,data['d2'][mask]+data['d1'][mask],'g',linewidth=0.6,label='d1+d2')
# ax_dist[0].legend()
# ax_dist[0].set_xlabel('Time')
# ax_dist[0].set_ylabel('nm')
# ax_dist[0].set_title("MetaD CVs")
#
# ax_dist[1].plot(timesteps_array,data['Min dist'][mask],'k',linewidth=0.6,label='min distance')
# ax_dist[1].legend()
# ax_dist[1].set_ylabel('nm')
# ax_dist[1].set_title("Minimum distance between the two chains")
# fig_dist.suptitle('Distance measurements')
# plt.show()
# sys.exit()


#Plotting CVs d1 and d2
fig_dist,ax_dist = plt.subplots(2,1,sharex=True)
ax_dist[0].plot(timesteps_array,data['d1'],'r',linewidth=2.6,label='d1')
ax_dist[0].plot(timesteps_array,data['d2'],'b',linewidth=2.6,label='d2')
#ax_dist[0].plot(timesteps_array,data['d2']+data['d1'],'g',linewidth=0.6,label='d1+d2')
ax_dist[0].legend(fontsize=20)
ax_dist[0].set_xlabel('Time',fontsize=40)
ax_dist[0].set_ylabel('nm',fontsize=40)
#ax_dist[0].set_title("MetaD CVs")
ax_dist[0].tick_params(labelsize=30)
ax_dist[1].plot(timesteps_array,data['Min dist'],'k',linewidth=0.6,label='min distance')
ax_dist[1].legend(fontsize=30)
ax_dist[1].set_ylabel('nm')
#ax_dist[1].set_title("Minimum distance between the two chains")
fig_dist.suptitle('Distance measurements')
plt.show()

#Plotting angles
fig_ang,ax_ang = plt.subplots()

ax_ang.plot(timesteps_array,data['Ax-ang'],'orange',linewidth=0.6,label='Axial Angle')
ax_ang.legend()
ax_ang.set_ylabel('Degrees')
ax_ang.set_title("Angle between long axes")
fig_ang.suptitle("Angle measurements")

#Plotting RMSD
mask=timesteps_array>0
fid_rmsd,ax_rmsd = plt.subplots(figsize=(12,8))
ax_rmsd.plot(timesteps_array[mask],data['RMSD-A'][mask],'tomato',linewidth=2.0,label='chain A')
ax_rmsd.plot(timesteps_array[mask],data['RMSD-BtoA'][mask],'olive',linewidth=2.0,label='chain B')
ax_rmsd.axhline(y=2.0,linestyle='--',linewidth=3.2)
ax_rmsd.legend(fontsize=30)
ax_rmsd.set_xlabel(r'Time ($ \mu s$)',fontsize=40)
ax_rmsd.set_ylabel('RMSD (nm)',fontsize=40)
#ax_rmsd.set_title("RMSD")
ax_rmsd.tick_params(labelsize=35)
# ax_rmsd.ticklabel_format(style='sci',fontsize=30)
fid_rmsd.tight_layout()
plt.show()


#Plotting d1 and d2
mask=timesteps_array>0
fid_rmsd,ax_rmsd = plt.subplots(figsize=(12,8))
ax_rmsd.plot(timesteps_array[mask],data['d1'][mask],'red',linewidth=2.0,label='d1')
ax_rmsd.plot(timesteps_array[mask],data['d2'][mask],'blue',linewidth=2.0,label='d2')
# ax_rmsd.axhline(y=2.0,linestyle='--',linewidth=3.2)
ax_rmsd.legend(fontsize=30)
ax_rmsd.set_xlabel(r'Time ($ \mu s$)',fontsize=40)
ax_rmsd.set_ylabel('CV distance (nm)',fontsize=40)
#ax_rmsd.set_title("RMSD")
ax_rmsd.tick_params(labelsize=35)
# ax_rmsd.ticklabel_format(style='sci',fontsize=30)
fid_rmsd.tight_layout()
plt.show()
