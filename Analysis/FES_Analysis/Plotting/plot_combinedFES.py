import matplotlib
import seaborn as sns
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import numpy as np
import sys
from matplotlib.contour import QuadContourSet as cs
# Get data

print(matplotlib.__version__)


def create_fesgrid(datafile):
    print("Reading data...")
    data = [line.split() for line in open(datafile, "r")]
    data2 = [x for x in data if not x == []]  # strip out headers
    file_name = str(sys.argv[1])
    d1, d2, free, dd1, dd2 = [], [], [], [], []
    for elem in data2[9:]:
        if elem[2] == 'inf' or elem[2] == '-inf':
            free.append(-10.0)
        else:
            free.append(float(elem[2]))

        d1.append(float(elem[0]))
        d2.append(float(elem[1]))
    #    dd1.append(float(elem[3]))
    #    dd2.append(float(elem[4]))

    X = np.linspace(min(d1), max(d1), 1318)
    Y = np.linspace(min(d2), max(d2), 1322)

    print("Creating data grid. This may take a while...")
    D1, D2 = np.meshgrid(X, Y)

    #Normalize unbound state to zero
    free = np.array(free)
    # free = free+110.6

    #Zero value
    d1_arr = np.array(d1)
    d2_arr = np.array(d2)
    mask1 = (d1_arr >= 15.0 ) & (d1_arr < 16.0 )
    mask2 = (d2_arr >= 15.0 ) & (d2_arr < 16.0 )
    mask3 = mask1 * mask2

    correction = np.mean(free[mask3])
    free = free - correction

    ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)

    #Shift max value to a constant. To shifht the color map
    # max_val = 5
    # mask4 = (free>=max_val)
    # free[mask4]=max_val

    return(D1,D2, ENER, free)


font_size=25
cbar_font_size=15
t_size=16
leg_size=20


fig,[ax,ax1,ax2]=plt.subplots(3,1,sharex=True,figsize=(6,34))

"""
Solution FES
"""
sol_file = sys.argv[1]
D1_sol,D2_sol, ENER_sol, free_sol = create_fesgrid(sol_file)

max_val=45
levels_sol = np.arange(np.min(free_sol),max_val, (max_val-np.min(free_sol))/30)
levels_sol = levels_sol.tolist()
levels_sol = levels_sol
levels_sol = np.array(levels_sol)

clr_map = cm.jet_r
clr_map.set_over('black')
max_color=20
contour = ax.contour(D1_sol, D2_sol, ENER_sol, colors='k', linewidths=0.3, levels=levels_sol)
contourf = ax.contourf(D1_sol, D2_sol, ENER_sol, cmap=clr_map,levels=levels_sol,alpha=0.9,vmin=np.min(free_sol),vmax=max_color)

# new_levels = contourf.levels[contourf.levels<60]
# contourf.set_clim(np.min(free_sol),60)

cbar = plt.colorbar(contourf,ax=ax,format="%4d", ticks=[np.min(free_sol),-25,0,25,max_val])
ax.scatter([0.8525], [0.7766], marker='x', c='forestgreen', s=300,alpha=1.0)
# ax.set_xlabel(r"$d_1$ (nm)",fontdict={'FontSize':font_size})
ax.set_ylabel(r"$d_2$ (nm)",fontdict={'fontsize':font_size})

cbar.set_label(label="Free Energy (kJ/mol)",size=cbar_font_size,weight='bold')
ax.tick_params(axis='y',labelsize=t_size)
cbar.ax.tick_params(axis='y',labelsize=t_size)
ax.set_xlim(left=0.5,right=19.0)
ax.set_ylim(bottom=0.5,top=19.0)


"""
Membrane FES
"""
mem_file = sys.argv[2]
D1_mem,D2_mem, ENER_mem, free_mem = create_fesgrid(mem_file)

levels_mem = np.arange(np.min(free_mem),np.max(free_mem), (np.max(free_mem)-np.min(free_mem))/30)
levels_mem = levels_mem.tolist()
levels_mem = levels_mem
levels_mem = np.array(levels_mem)

clr_map = cm.jet_r
clr_map.set_over('black')
max_color=25
contour = ax1.contour(D1_mem, D2_mem, ENER_mem, colors='k', linewidths=0.3, levels=levels_mem)
contourf = ax1.contourf(D1_mem, D2_mem, ENER_mem, cmap=clr_map,levels=levels_mem,alpha=0.9,vmax=max_color)


cbar = plt.colorbar(contourf,ax=ax1,format="%4d", ticks=[np.min(free_mem),-100,-75,-50,-25,0,25,50,75,100])
ax1.scatter([0.8525], [0.7766], marker='x', c='forestgreen', s=300,alpha=1.0)

ax1.set_ylabel(r"$d_2$ (nm)",fontdict={'fontsize':font_size})

cbar.set_label(label="Free Energy (kJ/mol)",size=cbar_font_size,weight='bold')
ax1.tick_params(axis='y',labelsize=t_size)
cbar.ax.tick_params(axis='y',labelsize=t_size)
ax1.set_xlim(left=0.5,right=19.0)
ax1.set_ylim(bottom=0.5,top=19.0)


"""
Pseudo FES
"""
ps_file = sys.argv[3]
D1_ps,D2_ps, ENER_ps, free_ps = create_fesgrid(ps_file)

levels_ps = np.arange(np.min(free_ps),np.max(free_ps), (np.max(free_ps)-np.min(free_ps))/30)
levels_ps = levels_ps.tolist()
levels_ps = levels_ps
levels_ps = np.array(levels_ps)

clr_map = cm.jet_r
clr_map.set_over('black')
max_color=25

contour = ax2.contour(D1_ps, D2_ps, ENER_ps, colors='k', linewidths=0.3, levels=levels_ps)
contourf = ax2.contourf(D1_ps, D2_ps, ENER_ps, cmap=clr_map,levels=levels_ps,alpha=0.9,vmax=max_color)

cbar = plt.colorbar(contourf,ax=ax2,format="%4d", ticks=[np.min(free_ps),-50,-25,0,25,50,75])
ax2.scatter([0.8525], [0.7766], marker='x', c='forestgreen', s=300,alpha=1.0)
ax2.set_xlabel(r"$d_1$ (nm)",fontdict={'fontsize':font_size})
ax2.set_ylabel(r"$d_2$ (nm)",fontdict={'fontsize':font_size})

cbar.set_label(label="Free Energy (kJ/mol)",size=cbar_font_size,weight='bold')
ax2.tick_params(axis='both',labelsize=t_size)
cbar.ax.tick_params(axis='y',labelsize=t_size)
ax2.set_xlim(left=0.5,right=19.0)
ax2.set_ylim(bottom=0.5,top=19.0)

ax2.set_xticks([5,10,15])

top=0.985
bottom=0.077
left=0.159
right=0.932
hspace=0.08
wspace=0.195
fig.tight_layout()
plt.subplots_adjust(top=top,bottom=bottom,left=left,right=right,hspace=hspace,wspace=wspace)
plt.savefig("FES_FINAL.svg",dpi=600,format="svg",transparent=True)
plt.savefig("FES_FINAL.png",dpi=600,format="png",transparent=True)

plt.show()


