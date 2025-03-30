import numpy as np 
import pandas
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.collections import LineCollection

file=sys.argv[1]
data = pandas.read_csv(file,comment='#',sep='\t')

timesteps_array = data['Timestep']/1e6


mask_bound = data['RMSD-BtoA'] < 2.0

lw = 2.5
f_size = 40
t_size = 36
msize = 17
fig_size=(10,8)


disp_ts = 100   #ns
interval = int(disp_ts*1e3 / 200)  #200 ps time interval between data points in the data file
filter_data = data[::interval] 


def get_segments(xdata,ydata):

    points = np.array([xdata, ydata]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

mask1 = (filter_data['d1'] < 12.0) & (filter_data['d2'] < 12.0)


"""
Plot d1
"""
color_map = np.where(mask1 , 'navy', 'silver')
line_maps = np.where(mask1 , 'solid', 'solid')
alpha_map = np.where(mask1 , 0.8, 0.8)
segments_d1 = get_segments(filter_data['Timestep']/1e6,filter_data['d1'])
lc_d1 = LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps)
fig,[ax,ax1]=plt.subplots(2,1,figsize=fig_size,sharex=True)
# ax.plot(data['Timestep'][::interval]/1e6,data['d1'][::interval],'crimson',linewidth=2.0,label=r'$d_1$',linestyle='-')
ax.add_collection(lc_d1)
# ax.legend(fontsize=30)
# ax.set_xlabel(r'Time ($ \mu s$)',fontsize=f_size)
ax.set_ylabel(r'$d_1$ (nm)',fontsize=f_size)
#ax_rmsd.set_title("RMSD")
ax.tick_params(labelsize=t_size)
ax.set_yticks([0,10,20,30])
# ax.yaxis.tick_right()
# ax_rmsd.ticklabel_format(style='sci',fontsize=30)


bd_data = data[::interval]
bd_data = bd_data[mask_bound]
ax.plot(bd_data['Timestep']/1e6,bd_data['d1'],color='black',linewidth=lw,label=r'$d_1$ bound',linestyle='',marker='*',markersize=msize)
# fig.tight_layout()
# figlbl = "d1_FINAL_" + str(disp_ts) 
# plt.savefig(figlbl+ "ns.svg",dpi=600,format="svg",transparent=True)
# plt.savefig(figlbl+ "ns.png",dpi=600,format="png",transparent=True)


"""
Plot d2
"""
color_map = np.where(mask1 , 'navy', 'silver')
line_maps = np.where(mask1 , 'solid', 'solid')
alpha_map = np.where(mask1 , 0.8, 0.8)
# fig,ax=plt.subplots(figsize=fig_size)
segments_d2 = get_segments(filter_data['Timestep']/1e6,filter_data['d2'])
lc_d2 = LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps)

# ax.plot(data['Timestep'][::interval]/1e6,data['d1'][::interval],'crimson',linewidth=2.0,label=r'$d_1$',linestyle='-')
ax1.add_collection(lc_d2)


ax1.plot(bd_data['Timestep']/1e6,bd_data['d2'],color='black',linewidth=lw,label=r'$d_2$ bound',linestyle='',marker='*',markersize=msize)
# ax.legend(fontsize=30)
ax1.set_xlabel(r'Time ($ \mu s$)',fontsize=f_size)
ax1.set_ylabel(r'$d_2$ (nm)',fontsize=f_size)
#ax_rmsd.set_title("RMSD")
ax1.tick_params(labelsize=t_size)
ax1.set_yticks([0,10,20,30])
# ax.yaxis.tick_right()
# ax_rmsd.ticklabel_format(style='sci',fontsize=30)
fig.tight_layout()
figlbl = "d1d2_FINAL_" + str(disp_ts) 
plt.savefig(figlbl+ "ns.svg",dpi=600,format="svg",transparent=True)
plt.savefig(figlbl+ "ns.png",dpi=600,format="png",transparent=True)

plt.show()


"""
Plot RMSD
"""
# fig,ax = plt.subplots(figsize=(12,8))

# # ax.plot(data['Timestep'][::interval]/1e6,data['d1'][::interval],'crimson',linewidth=2.0,label=r'$d_1$',linestyle='-')
# ax.plot(filter_data['Timestep']/1e6,filter_data['RMSD-B'],linewidth=2.0,linestyle='-',color='darkcyan',alpha=0.8)

# # ax.plot(bd_data['Timestep']/1e6,bd_data['d2'],color='black',linewidth=2.0,linestyle='',marker='*',markersize=10)
# # ax.legend(fontsize=30)
# ax.set_xlabel(r'Time ($ \mu s$)',fontsize=40)
# ax.set_ylabel(r'$RMSD$ (nm)',fontsize=40)
# ax.axhline(y=2.0,linestyle='-.',linewidth=2.0,color='black')
# ax.tick_params(labelsize=35)
# fig.tight_layout()
# figlbl = "RMSD_FINAL_" + str(disp_ts) + "ns.svg"
# plt.savefig(figlbl,dpi=600,format="svg",transparent=True)
# plt.show()

"""
Plot RMSD
"""
fig,ax = plt.subplots(figsize=fig_size)

# ax.plot(data['Timestep'][::interval]/1e6,data['d1'][::interval],'crimson',linewidth=2.0,label=r'$d_1$',linestyle='-')
ax.plot(filter_data['Timestep']/1e6,filter_data['comZ_chainA'],linewidth=lw,linestyle='-',color='sienna',alpha=0.8)

# ax.plot(bd_data['Timestep']/1e6,bd_data['d2'],color='black',linewidth=2.0,linestyle='',marker='*',markersize=10)
# ax.legend(fontsize=30)
ax.set_xlabel(r'Time ($ \mu s$)',fontsize=f_size)
ax.set_ylabel(r'$COM_{z} (nm)$',fontsize=f_size)
ax.set_ylim(0,36.0)

ax.set_yticks([0,9,18,27,36])
# ax.axhline(y=2.0,linestyle='-.',linewidth=2.0,color='black')
ax.tick_params(labelsize=t_size)
fig.tight_layout()
figlbl = "COM_FINAL_" + str(disp_ts) 
plt.savefig(figlbl+ "ns.svg",dpi=600,format="svg",transparent=True)
plt.savefig(figlbl+ "ns.png",dpi=600,format="png",transparent=True)
plt.show()

