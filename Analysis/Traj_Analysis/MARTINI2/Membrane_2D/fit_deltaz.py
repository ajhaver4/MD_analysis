import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm

data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")

n_frames = len(data['Timestep'])
timesteps_array = data['Timestep']/1e6
# timesteps_array = data['Frame']
threshold = data['Timestep'] > 0.1e5

fig,ax = plt.subplots()
ax.plot(timesteps_array,data['comA'],label='comA',linewidth=0.8,color='crimson')
ax.plot(timesteps_array,data['comB'],label='comB',linewidth=0.8,color='royalblue')
ax.set_ylabel("Distance (nm)")
ax.set_xlabel("Frame no. ")
ax.legend()

fig1,[ax1,ax2]=plt.subplots(2,1)
dzA=ax1.plot(timesteps_array,data['delta_z_A'],label='chain A',linewidth=0.8,color='crimson')
dzB=ax1.plot(timesteps_array,data['delta_z_B'],label='chain B',linewidth=0.8,color='royalblue')
fig1.suptitle("Delta Z values")
ax2.hist(data['delta_z_A'],bins=100,color='crimson',alpha=0.6,density=True)
ax2.hist(data['delta_z_B'],bins=100,color='royalblue',alpha=0.6,density=True)
ax1.set_ylabel("Distance (nm)")
ax1.set_xlabel("Frame no. ")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Distance (nm)")
ax1.legend()
ax2.legend()

fig3,[ax3,ax4,ax5]=plt.subplots(3,1)
ax3.hist(data['z_A_hi'],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax3.hist(data['z_B_hi'],bins=100,label='chain B',alpha=0.6,density=True,color='royalblue')
ax3.set_ylabel("Frequency")
ax3.set_xlabel("Distance (nm)")
ax3.set_title("High flexibility region")
ax3.legend()

ax4.hist(data['z_A_low'],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax4.hist(data['z_B_low'],bins=100,label='chain B',alpha=0.6,density=True,color='royalblue')
ax4.set_ylabel("Frequency")
ax4.set_xlabel("Distance (nm)")
ax4.set_title("Low flexibility region")
ax4.legend()

ax5.hist(data['z_A_med'],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax5.hist(data['z_B_med'],bins=100,label='chain B',alpha=0.6,density=True,color='royalblue')
ax5.set_ylabel("Frequency")
ax5.set_xlabel("Distance (nm)")
ax5.set_title("Mix flexibility region")
ax5.legend()

plt.show()


def harm2(x,k,r0):
    c = (2*3.14*2.5/k)**0.5
    v = 0.5*k*(x-r0)**2 + 2.5*np.log(c)
    return(v)


"""
Fitting membrane Distances
"""
figm,[axm1,axm2]=plt.subplots(2,1)
axm1.plot(timesteps_array[threshold],data['z_A_membres'][threshold],label='chain A',linewidth=0.8,color='crimson')
axm1.plot(timesteps_array[threshold],data['z_B_membres'][threshold],label='chain B',linewidth=0.8,color='steelblue')
fig1.suptitle("Membrane residues Delta Z values")

maskm = threshold & (data['z_A_membres'] < 1.3 ) & (data['z_A_membres'] > 0.7)  
memb_dist_A = data['z_A_membres'][maskm]

# maskb = (data['z_B_membres'] != 'nan') & threshold & (data['z_B_membres'] > 0.5 )
memb_dist_B = data['z_B_membres'][threshold]


pd_memb_A = axm2.hist(memb_dist_A,bins=100,color='crimson',alpha=0.6,density=True,label='ChainA')
pd_memb_B = axm2.hist(memb_dist_B,bins=100,color='steelblue',alpha=0.6,density=True,label='ChainB')


axm1.set_ylabel("Distance (nm)")
axm1.set_xlabel("Frame no. ")
axm2.set_ylabel("Frequency")
axm2.set_xlabel("Distance (nm)")
axm1.legend()
axm2.legend()

nzero_maskA = np.nonzero(pd_memb_A[0])
xdata_m = pd_memb_A[1][nzero_maskA]
ydata_m = pd_memb_A[0][nzero_maskA]


nzero_maskB = np.nonzero(pd_memb_B[0])
xdata_mB = pd_memb_B[1][nzero_maskB]
ydata_mB = pd_memb_B[0][nzero_maskB]

#Confidence Interval 
# Compute 95% confidence intervals
alpha = 0.05  # for 95% confidence
z_value = norm.ppf(1 - alpha / 2)  # Z value for 95% confidence level


params_m = curve_fit(harm2,xdata_m,-2.5*np.log(ydata_m),[334,1.25],bounds=((330,1.15),(338,1.35)))
params_mB = curve_fit(harm2,xdata_mB,-2.5*np.log(ydata_mB),[490,1.2],bounds=((488,1.15),(496,1.35)))

"Plotting parameters   "
t_size=22
f_size=23
m_size=12
leg_size=24
fig_size=(6,5)


print("Chain A dz fitting: ")
print("k = {:.2f} kJ mol-1 nm-2".format(params_m[0][0]))
print("mu = {:.2f} nm".format(params_m[0][1]))
print("Cov: ",params_m[1])
print("STD of error : ", np.sqrt(np.diag(params_m[1])))
ci = []
std_err_02A = np.sqrt(np.diag(params_m[1]))
for param, err in zip(params_m[0], std_err_02A):
    ci.append([param - z_value * err, param + z_value * err])
print("95% CI : ", ci)


print("Chain B dz fitting: ")
print("k = {:.2f} kJ mol-1 nm-2".format(params_mB[0][0]))
print("mu = {:.2f} nm".format(params_mB[0][1]))
print("Cov: ",params_mB[1])
print("STD of error : ", np.sqrt(np.diag(params_mB[1])))
ci = []
std_err_02B = np.sqrt(np.diag(params_mB[1]))
for param, err in zip(params_mB[0], std_err_02B):
    ci.append([param - z_value * err, param + z_value * err])
print("95% CI : ", ci)

figm3,axm3 = plt.subplots(figsize=fig_size)
axm3.plot(xdata_m,-2.5*np.log(ydata_m),label=r'$U_{2D} (dz)$',alpha=0.7,color='skyblue',linestyle='dashed',marker='o',markersize=m_size)
axm3.plot(xdata_m,harm2(xdata_m,params_m[0][0],params_m[0][1]),label=r'$U_{fit} (dz)$',alpha=0.7,color='k',linewidth=2.0)
axm3.legend(fontsize=leg_size,frameon=False)
# axm3.set_title("Chain A membrane Res")
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_m[0][0])
s2 = r'$mu$ = {:.2f}'.format(params_m[0][1])
# axm3.text(np.mean(data['z_A_membres']),10,s1)
# axm3.text(np.mean(data['z_A_membres']),11,s2)
figm3.tight_layout()

axm3.tick_params(labelsize=t_size)
axm3.set_xlabel(r'$z_{surf}$ (nm)',fontsize=f_size)
axm3.set_ylabel("Potential Energy (kJ/mol)",fontsize=f_size)


figm4,axm4 = plt.subplots(figsize=fig_size)
axm4.plot(xdata_mB,-2.5*np.log(ydata_mB),label=r'$U_{2D} (dz)$',alpha=0.7,color='skyblue',linestyle='dashed',marker='o',markersize=m_size)
axm4.plot(xdata_mB,harm2(xdata_mB,params_mB[0][0],params_mB[0][1]),label=r'$U_{fit} (dz)$',alpha=0.7,color='k',linewidth=2.0)
axm4.legend(fontsize=leg_size,frameon=False)
# axm4.set_title("Chain B membrane Res")
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_mB[0][0])
s2 = r'$mu$ = {:.2f}'.format(params_mB[0][1])
# axm4.text(np.mean(data['z_B_membres']),10,s1)
# axm4.text(np.mean(data['z_B_membres']),11,s2)
figm4.tight_layout()
axm4.tick_params(labelsize=t_size)
# axm4.set_ylim(0,4.5)
axm4.set_xlabel(r'$z_{surf}$ (nm)',fontsize=f_size)
axm4.set_ylabel("Potential Energy (kJ/mol)",fontsize=f_size)



k = np.linspace(10,200,10)
x=np.linspace(-2,2,50)


fig,ax = plt.subplots()
for i in range(len(k)):
    u=0.5*k[i]*(x**2)
    ax.plot(x,u,label=k[i],alpha=0.7,linewidth=0.8)

ax.legend()
ax.set_xlabel('r (nm)')
ax.set_ylabel('U(r)')
plt.show()


fig,ax = plt.subplots(figsize=fig_size)

pd_memb_A = ax.hist(memb_dist_A,bins=100,label='chain A',alpha=0.6,density=True,color='skyblue')
fig.tight_layout()
ax.tick_params(labelsize=t_size)
ax.set_ylim(0,6)
ax.set_xlabel(r'$z_{surf}$ (nm)',fontsize=f_size)
ax.set_ylabel(r"$P_{2D} (z_{surf})$",fontsize=f_size)

fig,ax = plt.subplots(figsize=fig_size)
pd_memb_B = ax.hist(memb_dist_B,bins=100,color='skyblue',alpha=0.6,density=True,label='ChainB')
fig.tight_layout()
ax.tick_params(labelsize=t_size)
ax.set_ylim(0,6)
ax.set_xlabel(r'$z_{surf}$ (nm)',fontsize=f_size)
ax.set_ylabel(r"$P_{2D} (z_{surf})$",fontsize=f_size)

plt.show()