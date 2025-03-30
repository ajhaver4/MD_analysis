import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm


"Plotting parameters   "
t_size=22
f_size=23
m_size=12
leg_size=24
fig_size=(6,5)

# col = 'Frame\tTimestep\tcomA\tdelta_z_A\tz_A_hi\tz_A_low\tz_A_med\tz_A_membres'
# col_names=col.split("\t")
# print("Column Names: ",col_names)
# data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t",header=None,names=col_names)
z_data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")
final_data = pd.read_csv(sys.argv[2],comment='#',sep='\t')

#BOUND STATE FILTERING
filter_time = []
for i in range(len(final_data['d1'])):
    if final_data['d1'][i]<2.0 or final_data['d1'][i]<2.0:
        filter_time.append(final_data['Timestep'][i])
data=z_data[z_data['Timestep'].isin(filter_time)]


#FOR UNBOUND STATE FILTERING
# filter_time_A=[]
# filter_time_B=[]
# for i in range(len(final_data['d1'])):
#     if final_data['d1'][i]>12.0 and final_data['d1'][i]>12.0:
#         if final_data['min_dist_A'][i] <= 1.0:
#             filter_time_A.append(final_data['Timestep'][i])
#         if final_data['min_dist_B'][i] <=1.0:
#             filter_time_B.append(final_data['Timestep'][i])

dataA = z_data[z_data['Timestep'].isin(filter_time)]
dataB = z_data[z_data['Timestep'].isin(filter_time)]

# dataA = z_data
# dataB = z_data
n_frames = len(dataA['Timestep'])
timesteps_arrayA = dataA['Timestep']/1e6
timesteps_arrayB = dataB['Timestep']/1e6


# thresholdA = (dataA['Timestep'] > 0.6*1e6) | (dataA['Timestep'] < 0.2*1e6)
# thresholdB = (dataB['Timestep'] > 0.6*1e6) | (dataB['Timestep'] < 0.2*1e6)

thresholdA = dataA['Timestep'] > 0.2e6
thresholdB = dataB['Timestep'] > 0.2e6

# print(thresholdA)
# sys.exit()
# timesteps_array = data['Frame']

fig,ax = plt.subplots()
ax.plot(timesteps_arrayA,dataA['comA'],label='comA',linewidth=0.8,color='crimson')
ax.set_ylabel("Distance (nm)")
ax.set_xlabel("Frame no. ")
ax.legend()

fig1,[ax1,ax2]=plt.subplots(2,1)
ax1.plot(timesteps_arrayA,dataA['delta_z_A'],label='chain A',linewidth=0.8,color='crimson')
fig1.suptitle("Delta Z values")
pd_zA = ax2.hist(dataA['delta_z_A'][thresholdA],bins=100,color='crimson',alpha=0.6,density=True,label='ChainA')
ax1.set_ylabel("Distance (nm)")
ax1.set_xlabel("Frame no. ")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Distance (nm)")
ax1.legend()
ax2.legend()


fig2,[ax1,ax2]=plt.subplots(2,1)
ax1.plot(timesteps_arrayB,dataB['delta_z_B'],label='chain B',linewidth=0.8,color='steelblue')
fig1.suptitle("Delta Z values")
pd_zB = ax2.hist(dataB['delta_z_B'][thresholdB],bins=100,color='steelblue',alpha=0.6,density=True,label='ChainB')
ax1.set_ylabel("Distance (nm)")
ax1.set_xlabel("Frame no. ")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Distance (nm)")
ax1.legend()



plt.show()


#Curve fitting
def harm(x,k,r0):
    v = 0.5*k*(x-r0)**2
    return(v)

def harm2(x,k,r0):
    c = (2*3.14*2.5/k)**0.5
    v = 0.5*k*(x-r0)**2 + 2.5*np.log(c)
    return(v)

# print(len(pd_hi_A[1]))
# print(len(pd_hi_A[0]))
# sys.exit()


print("Statistics for Chain A: ")
mean_a = np.mean(dataA['delta_z_A'])
std_a = np.std(dataA['delta_z_A'])
print("Mean: ",mean_a)
print("STD: ",std_a)

print("Statistics for Chain B: ")
mean_b = np.mean(dataB['delta_z_B'])
std_b = np.std(dataB['delta_z_B'])
print("Mean: ",mean_b)
print("STD: ",std_b)



mask1 = np.nonzero(pd_zA[0])
mask2 = np.nonzero(pd_zB[0])

params_A = curve_fit(harm2,pd_zA[1][mask1],-2.5*np.log(pd_zA[0][mask1]),[12,1.0])
params_B = curve_fit(harm2,pd_zB[1][mask2],-2.5*np.log(pd_zB[0][mask2]),[12,1.2])


#Confidence Interval 
# Compute 95% confidence intervals
alpha = 0.05  # for 95% confidence
z_value = norm.ppf(1 - alpha / 2)  # Z value for 95% confidence level

print("Chain A COM dz fitting: ")
print("k = {:.2f} kJ mol-1 nm-2".format(params_A[0][0]))
print("mu = {:.2f} nm".format(params_A[0][1]))
print("Cov: ",params_A[1])
print("STD of error : ", np.sqrt(np.diag(params_A[1])))
ci = []
std_err_02A = np.sqrt(np.diag(params_A[1]))
for param, err in zip(params_A[0], std_err_02A):
    ci.append([param - z_value * err, param + z_value * err])
print("95% CI : ", ci)


print("Chain B COM dz fitting: ")
print("k = {:.2f} kJ mol-1 nm-2".format(params_B[0][0]))
print("mu = {:.2f} nm".format(params_B[0][1]))
print("Cov: ",params_B[1])
print("STD of error : ", np.sqrt(np.diag(params_B[1])))
ci = []
std_err_02B = np.sqrt(np.diag(params_B[1]))
for param, err in zip(params_B[0], std_err_02B):
    ci.append([param - z_value * err, param + z_value * err])
print("95% CI : ", ci)


fig_fit,[ax_fit2,ax_fit3] = plt.subplots(2,1)

ax_fit2.plot(pd_zA[1][mask1],-2.5*np.log(pd_zA[0][mask1]),label='V(x)=-2.5 * log(P(x))')
ax_fit2.plot(pd_zA[1][mask1],harm2(pd_zA[1][mask1],*params_A[0]),label='V(x) = 0.5*k*x^2')
ax_fit2.legend()
ax_fit2.set_title('DeltaZ for Chain A')
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_A[0][0])
ax_fit2.text(np.mean(dataA['z_A_membres']),10,s1)

ax_fit3.plot(pd_zB[1][mask2],-2.5*np.log(pd_zB[0][mask2]),label='V(x)=-2.5 * log(P(x))')
ax_fit3.plot(pd_zB[1][mask2],harm2(pd_zB[1][mask2],*params_B[0]),label='V(x) = 0.5*k*x^2')
ax_fit3.legend()
ax_fit3.set_title('DeltaZ for Chain B')
fig_fit.suptitle("Chain B Fitting")
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_B[0][0])
ax_fit3.text(np.mean(dataB['z_B_membres']),10,s1)



plt.show()


"""
Fitting membrane Distances
"""
figm,[axm1,axm2]=plt.subplots(2,1)
axm1.plot(timesteps_arrayA[thresholdA],dataA['z_A_membres'][thresholdA],label='chain A',linewidth=0.8,color='crimson')
axm1.plot(timesteps_arrayB[thresholdB],dataA['z_B_membres'][thresholdB],label='chain B',linewidth=0.8,color='steelblue')
fig1.suptitle("Membrane residues Delta Z values")

maskm = (dataA['z_A_membres'] != 'nan') & thresholdA & (dataA['z_A_membres'] > 0.2 )
memb_dist_A = dataA['z_A_membres'][maskm]

# maskb = (dataA['z_B_membres'] != 'nan') & thresholdB & (dataA['z_B_membres'] > 0.5) & (dataA['z_B_membres'] < 1.8)
maskb = (dataA['z_B_membres'] != 'nan') & thresholdB & (dataA['z_B_membres'] > 0.8) & (dataA['z_B_membres'] < 1.9)
memb_dist_B = dataA['z_B_membres'][maskb]


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
# ydata_m = pd_memb_A[0][nzero_maskA]*(pd_memb_A[1][1]-pd_memb_A[1][0])

# ydata_mB = pd_memb_B[0][nzero_maskB]*(pd_memb_B[1][1]-pd_memb_B[1][0])
ydata_m = pd_memb_A[0][nzero_maskA]

# print(sum(pd_lo_A[0])*(pd_lo_A[1][1]-pd_lo_A[1][0]))
# print(ydata_m)
# sys.exit()
nzero_maskB = np.nonzero(pd_memb_B[0])
xdata_mB = pd_memb_B[1][nzero_maskB]
ydata_mB = pd_memb_B[0][nzero_maskB]
# ydata_mB = pd_memb_B[0][nzero_maskB]*(pd_memb_B[1][1]-pd_memb_B[1][0])

params_m = curve_fit(harm2,xdata_m,-2.5*np.log(ydata_m),[160,1.2])
params_mB = curve_fit(harm2,xdata_mB,-2.5*np.log(ydata_mB),[150,1.2],bounds=((140,1.15),(160,1.35)))


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
ax.set_ylim(0,4.5)
ax.set_xlabel(r'$z_{surf}$ (nm)',fontsize=f_size)
ax.set_ylabel(r"$P_{2D} (z_{surf})$",fontsize=f_size)

fig,ax = plt.subplots(figsize=fig_size)
pd_memb_B = ax.hist(memb_dist_B,bins=100,color='skyblue',alpha=0.6,density=True,label='ChainB')
fig.tight_layout()
ax.tick_params(labelsize=t_size)
ax.set_ylim(0,4.5)
ax.set_xlabel(r'$z_{surf}$ (nm)',fontsize=f_size)
ax.set_ylabel(r"$P_{2D} (z_{surf})$",fontsize=f_size)

plt.show()