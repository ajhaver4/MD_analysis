import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import constants
from scipy.stats import norm

# data = pd.read_csv(sys.argv[1],delimiter="\s+",names=['d1','fes'],comment='#')
data_d1 = pd.read_csv(sys.argv[1],delimiter="\s+",names=['d1','fes','std','err'],comment='#')
data_d2 = pd.read_csv(sys.argv[2],delimiter="\s+",names=['d1','fes','std','err'],comment='#')
data_d3 = pd.read_csv(sys.argv[3],delimiter="\s+",names=['d1','fes','std','err'],comment='#')
data_d4 = pd.read_csv(sys.argv[4],delimiter="\s+",names=['d1','fes','std','err'],comment='#')
data_d5 = pd.read_csv(sys.argv[5],delimiter="\s+",names=['d1','fes','std','err'],comment='#')
data_d6 = pd.read_csv(sys.argv[6],delimiter="\s+",names=['d1','fes','std','err'],comment='#')

cf = 0.95    # % CI
# norm_error = norm.ppf((cf+1)/2)
norm_error =1 

def plot_1Dfes(data,label,ax,alpha=1.0,color='royalblue',lw=1,ls='-'):

    mask_zero = (data['d1'] > 12.0) & (data['d1'] < 18.0)
    zero_val = np.mean(data['fes'][mask_zero])
    data['fes'] = data['fes'] - zero_val
    ax.plot(data['d1'],data['fes'],label=label,alpha=alpha,color=color,linewidth=lw,linestyle=ls)
    ax.fill_between(data['d1'],data['fes']+norm_error*data['err'],data['fes']-norm_error*data['err'],color=color,alpha=0.3)
    # ax.errorbar(data['d1'],data['fes'],yerr=data['err'],color=color,linewidth=lw,linestyle=ls,capsize=3,capthick=1)





#Order of inputs 
#3D_fes_d1.dat 3D_fes_d2.dat 2D_fes_d1.dat 2D_fes_d2.dat
input_labels=['3D_d1','3D_d2','2D_d1','2D_d2','PS_d1','PS_d2']


fig_size=(6,4)
f_size=25
t_size=20
leg_size = 20
lw=2.4
fig,ax1 = plt.subplots(figsize=fig_size)

ylim = -120



plot_1Dfes(data_d1,'3D',ax1,alpha=0.8,color='navy',lw=lw)
plot_1Dfes(data_d3,'2D',ax1,alpha=0.8,color='firebrick',lw=lw)
plot_1Dfes(data_d5,'2D-Ps',ax1,alpha=0.8,color='teal',lw=lw)

ax1.set_xlabel(r'$d_1 (nm)$',fontsize=f_size)    
ax1.set_ylabel(r"$F(d_1)$ (kJ/mol)",fontsize=f_size)
ax1.set_xlim(left=0.4,right=19.0)
ax1.set_ylim(top=25,bottom=ylim)
ax1.tick_params(labelsize=t_size)
ax1.legend(fontsize=leg_size,frameon=False)
fig.tight_layout()

ax1.vlines(0.81,np.min(data_d3['fes']),0,linestyle='--',color='black',linewidth=2.0)

plt.savefig("FES_d1_FINAL.svg",dpi=600,format="svg",transparent=True)
plt.savefig("FES_d1_FINAL.png",dpi=600,format="png",transparent=True)

fig,ax2 = plt.subplots(figsize=fig_size)
plot_1Dfes(data_d2,'3D',ax2,alpha=0.8,color='navy',lw=lw)
plot_1Dfes(data_d4,'2D',ax2,alpha=0.8,color='firebrick',lw=lw)
plot_1Dfes(data_d6,'2D-Ps',ax2,alpha=0.8,color='teal',lw=lw)

ax2.set_xlabel(r'$d_2 (nm)$',fontsize=f_size)    
ax2.set_ylabel(r"$F(d_2)$ (kJ/mol)",fontsize=f_size)
ax2.set_xlim(left=0.4,right=19.0)
ax2.set_ylim(top=25,bottom=ylim)
ax2.tick_params(labelsize=t_size)
# ax2.legend(fontsize=leg_size)
ax2.vlines(0.75,np.min(data_d4['fes']),0,linestyle='--',color='black',linewidth=2.0)
fig.tight_layout()
plt.savefig("FES_d2_FINAL.svg",dpi=600,format="svg",transparent=True)
plt.savefig("FES_d2_FINAL.png",dpi=600,format="png",transparent=True)

plt.show()


def integrate_fes(x_val,y_val,boundary=9.0):
    k = constants.value('Boltzmann constant')
    Ava_no = constants.value('Avogadro constant')
    temp = 310
    kbt = (k*temp*Ava_no)/1000   #kJ/mol
    fes = y_val
    rmsd = x_val
    
    dr = rmsd[1]-rmsd[0]

    mask_bd = (rmsd < boundary)
    mask_unb = (rmsd >= boundary)
    
    fes_bd_bins = fes[mask_bd].to_numpy()
    fes_unb_bins = fes[mask_unb].to_numpy()
 

    Z_bound = np.sum(np.exp(-fes_bd_bins/kbt)*dr)
    Z_unb = np.sum(np.exp(-fes_unb_bins/kbt)*dr)

    dG = -kbt*np.log(Z_bound/Z_unb)


    return dG


dG_3D_d1 = integrate_fes(data_d1['d1'],data_d1['fes'])
dG_3D_d2 = integrate_fes(data_d2['d1'],data_d2['fes'])
dG_2D_d1 = integrate_fes(data_d3['d1'],data_d3['fes'])
dG_2D_d2 = integrate_fes(data_d4['d1'],data_d4['fes'])


print("dG 3D d1: ",dG_3D_d1)
print("dG 3D d2: ",dG_3D_d2)
print("dG 2D d1: ",dG_2D_d1)
print("dG 2D d2: ",dG_2D_d2)


""""
Plotting both d1 and d2 FES on same plots for each environment
"""

#3D
fig_size=(6,10)
f_size=36
t_size=32
leg_size = 32
lw=4.2

fig,[ax1,ax2,ax3] = plt.subplots(3,1,figsize=fig_size,sharex=True)
plot_1Dfes(data_d1,r'$d_1$',ax1,alpha=0.8,color='navy',lw=lw)
plot_1Dfes(data_d2,r'$d_2$',ax1,alpha=0.8,color='navy',lw=lw,ls='dotted')

# ax1.set_xlabel(r'$\vec{s} (nm)$',fontsize=f_size)    
ax1.set_ylabel(r"$F_{3D}(\vec{s})$",fontsize=f_size)
ax1.set_xlim(left=0.4,right=19.0)
ax1.set_ylim(bottom=ylim,top=25)
ax1.tick_params(labelsize=t_size)
ax1.legend(fontsize=leg_size,loc ='lower right',frameon=False,labelspacing=0.2,handlelength=0.8)


ax1.vlines(0.81,-110,0,linestyle='--',color='black',linewidth=2.0)

# plt.savefig("FES_3D_only.svg",dpi=600,format="svg",transparent=True)
# plt.savefig("FES_3D_only.png",dpi=600,format="png",transparent=True)

#2D



plot_1Dfes(data_d3,r'$d_1$',ax2,alpha=0.8,color='firebrick',lw=lw)
plot_1Dfes(data_d4,r'$d_2$',ax2,alpha=0.8,color='firebrick',lw=lw,ls='dotted')

# ax2.set_xlabel(r'$\vec{s} (nm)$',fontsize=f_size)    
ax2.set_ylabel(r"$F_{2D}(\vec{s})$",fontsize=f_size)
ax2.set_xlim(left=0.4,right=19.0)
ax2.set_ylim(bottom=ylim,top=25)
ax2.tick_params(labelsize=t_size)
ax2.legend(fontsize=leg_size,loc ='lower right',frameon=False,labelspacing=0.2,handlelength=0.8)


ax2.vlines(0.81,-110,0,linestyle='--',color='black',linewidth=2.0)

# plt.savefig("FES_2D_only.svg",dpi=600,format="svg",transparent=True)
# plt.savefig("FES_2D_only.png",dpi=600,format="png",transparent=True)

#2D-Ps



plot_1Dfes(data_d5,r'$d_1$',ax3,alpha=0.8,color='teal',lw=lw)
plot_1Dfes(data_d6,r'$d_2$',ax3,alpha=0.8,color='teal',lw=lw,ls='dotted')

ax3.set_xlabel(r'$\vec{s} (nm)$',fontsize=f_size)    
ax3.set_ylabel(r"$F_{2DPs}(\vec{s})$",fontsize=f_size)
ax3.set_xlim(left=0.4,right=19.0)
ax3.set_ylim(bottom=ylim,top=25)

ax3.tick_params(labelsize=t_size)
ax3.legend(fontsize=leg_size,loc ='lower right',frameon=False,labelspacing=0.2,handlelength=0.8)


ax3.vlines(0.81,-110,0,linestyle='--',color='black',linewidth=2.0)
fig.tight_layout()
# plt.show()
plt.savefig("FES_Combined.svg",dpi=600,format="svg",transparent=True)
plt.savefig("FES_Combined.png",dpi=600,format="png",transparent=True)


#Output minimum energy values for each environment
def get_min_fes(data):
    mask_zero = (data['d1'] > 12.0) & (data['d1'] < 18.0)
    zero_val = np.mean(data['fes'][mask_zero])
    data['fes'] = data['fes'] - zero_val
    min_fes_indx = np.argmin(data['fes'].to_numpy())
    min_fes_d1 = data['d1'].to_numpy()[min_fes_indx]
    
    return (data['fes'][min_fes_indx],min_fes_d1)


dG_3D_d1 = get_min_fes(data_d1)
dG_3D_d2 = get_min_fes(data_d2)
dG_2D_d1 = get_min_fes(data_d3)
dG_2D_d2 = get_min_fes(data_d4)
dG_ps_d1 = get_min_fes(data_d5)
dG_ps_d2 = get_min_fes(data_d6)
print("dG 3D d1: ",dG_3D_d1)
print("dG 3D d2: ",dG_3D_d2)
print("dG 2D d1: ",dG_2D_d1)
print("dG 2D d2: ",dG_2D_d2)
print("dG ps d1: ",dG_ps_d1)
print("dG ps d2: ",dG_ps_d2)


