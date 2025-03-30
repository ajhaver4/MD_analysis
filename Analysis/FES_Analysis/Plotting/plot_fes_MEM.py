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

print("Reading data...")
data = [line.split() for line in open(sys.argv[1], "r")]
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
print(np.max(free))
print(np.min(free))
correction = np.mean(free[mask3])
free = free - correction

#Shift max value to a constant. To shifht the color map
# max_val = 5
# mask4 = (free>=max_val)
# free[mask4]=max_val


font_size=25
t_size=20
leg_size=20


fig=plt.figure(figsize=(8,6))
ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)

levels = np.arange(np.min(free),np.max(free), (np.max(free)-np.min(free))/30)
levels = levels.tolist()
levels = levels
levels = np.array(levels)
# print(levels)
# print(len(levels))
clr_map = cm.jet_r

# clr_map.set_under('white')
# clr_map.set_extremes('magenta')
contour = plt.contour(D1, D2, ENER, colors='k', linewidths=0.3, levels=levels)
contourf = plt.contourf(D1, D2, ENER, cmap=clr_map,levels=levels,alpha=0.9)
# print(contourf.levels)

# new_levels = contourf.levels[contourf.levels<60]
# print(new_levels)
# contourf.set_clim(np.min(free),60)
# con_set = cs(plt.gca(),D1, D2, ENER,levels= new_levels,cmap=clr_map)

cbar = plt.colorbar(contourf,format="%4d", ticks=[np.min(free),-100,-75,-50,-25,0,25,50,75])
plt.scatter([0.8525], [0.7766], marker='x', c='magenta', s=200)
plt.xlabel(r"$d_1$ (nm)",fontdict={'FontSize':font_size})
plt.ylabel(r"$d_2$ (nm)",fontdict={'FontSize':font_size})

ax=plt.gca()
cbar.set_label(label="Free Energy (kJ/mol)",size=font_size,weight='bold')
ax.tick_params(axis='both',labelsize=t_size)
cbar.ax.tick_params(axis='y',labelsize=t_size)
# cbar.set_clim(np.min(free),60)
# sns.kdeplot(d1,d2,cmap='Spectral',shade=True)
ax.set_yticks([0,5,10,15])
ax.set_xticks([0,5,10,15])
ax.set_xlim(left=0.5,right=19.0)
ax.set_ylim(bottom=0.5,top=19.0)
fig.tight_layout()
plt.savefig("FES_FINAL.svg",dpi=600,format="svg",transparent=True)
plt.savefig("FES_FINAL.png", dpi=900,format='png', bbox_inches='tight')
plt.show()
