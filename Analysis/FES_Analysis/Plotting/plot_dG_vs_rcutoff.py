import pandas
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

data3D = pandas.read_csv(sys.argv[1],comment='#',delimiter='\t')
data3D_m2 = pandas.read_csv(sys.argv[2],comment='#',delimiter='\t')


#3D, 3D_m2 and 3D_m2-Ps
avg_start_times = [104,60]
# avg_dG_values = [-23.341,-95.972,-49.434]

times3D = data3D['Timestep'].to_numpy()
times3D_m2 = data3D_m2['Timestep'].to_numpy()

indx3D = np.argmin(np.abs(times3D - avg_start_times[0]*1e6))
indx3D_m2 = np.argmin(np.abs(times3D_m2 - avg_start_times[1]*1e6))

avg_data3D = data3D[indx3D:]
avg_data3D_m2 = data3D_m2[indx3D_m2:]



ts_3D = times3D[1]-times3D[0]
ts_3D_m2 = times3D_m2[1]-times3D_m2[0]

print("Timestep 3D: ", ts_3D)
print("Timestep 3D_m2: ", ts_3D_m2)

block_size = 20   #us

time_int3D = int(block_size*1e6/ts_3D)
time_int3D_m2 = int(block_size*1e6/ts_3D_m2)

nblocks3D = int(len(avg_data3D['TS2'][::time_int3D]))
nblocks3D_m2 = int(len(avg_data3D_m2['TS2'][::time_int3D_m2]))

print("Number of blocks 3D: ", nblocks3D)
print("Number of blocks 3D_m2: ", nblocks3D_m2)
#Get the dG's for each time interval
# dG_blocks3D = {std_timestamps[i] : avg_data3D['TS2'][::time_int3D[i]] for i in range(len(std_timestamps))}
# dG_blocks3D_m2 = {std_timestamps[i] : avg_data3D_m2['TS2'][::time_int3D_m2[i]] for i in range(len(std_timestamps))}
# dG_blocksPs = {std_timestamps[i] : avg_dataPs['TS2'][::time_intPs[i]] for i in range(len(std_timestamps))}

std_values3D = {state: np.std(avg_data3D[state][::time_int3D],ddof=1) for state in avg_data3D.columns[1:]}
std_values3D_m2 = {state: np.std(avg_data3D_m2[state][::time_int3D_m2],ddof=1) for state in avg_data3D_m2.columns[1:]}



avg_values3D = {state : np.mean(avg_data3D[state][::time_int3D]) for state in avg_data3D.columns[1:]}
avg_values3D_m2 = {state : np.mean(avg_data3D_m2[state][::time_int3D_m2]) for state in avg_data3D_m2.columns[1:]}


cf = 0.95    # % CI
norm_error = norm.ppf((cf+1)/2)
error_bars3D = { state : norm_error*std_values3D[state]/np.sqrt(nblocks3D) for state in std_values3D.keys()}
error_bars3D_m2 = { state : norm_error*std_values3D_m2[state]/np.sqrt(nblocks3D_m2) for state in std_values3D_m2.keys()}




#For the plot line select which block size we wnat to use to plot the corresponding avg dG 
avg_dGs = [avg_values3D['TS2'],avg_values3D_m2['TS2']]



std_data3D = np.std(avg_data3D['TS2'])
std_data3D_m2 = np.std(avg_data3D_m2['TS2'])

print("Standard Deviation 3D: ", std_data3D)
print("Standard Deviation 3D_m2: ", std_data3D_m2)

print("--------------------------------------------------------")
print("Average Values 3D: ", avg_values3D)
print("Average Values 3D_m2: ", avg_values3D_m2)

print("----------------------------------------------------------")
print("Standard Deviation for a range of different time intervals: ")
print("Standard Deviation 3D: ", std_values3D)
print("Standard Deviation 3D_m2: ", std_values3D_m2)

print("-------------------------------------------------------------")
print("Error Bars for a range of different time intervals: ")
print("Error Bars 3D: ", error_bars3D)
print("Error Bars 3D_m2: ", error_bars3D_m2)

fig_size=(8,4)
f_size=20
t_size=15
leg_size = 14
lw=2.0
fig,[ax,ax1] = plt.subplots(1,2,figsize=fig_size)


x_values = [2.0,4.0,6.0,8.0,10.0,11.0,12.0,12.5,13.0,13.5,14.0,14.5,15.0,16.0]

ax.plot(x_values,avg_values3D.values(),linewidth=lw,alpha=0.8,label=r'$\Delta G_{3D} $',color='black')
ax1.plot(x_values,avg_values3D_m2.values(),linewidth=lw,alpha=0.8,label=r'$\Delta G_{3D} $',color='black')

ax.errorbar(x_values,avg_values3D.values(),yerr=error_bars3D.values(),fmt='o',color='black',capsize=5)
ax1.errorbar(x_values,avg_values3D_m2.values(),yerr=error_bars3D.values(),fmt='o',color='black',capsize=5)



ax.set_ylabel(r'$\Delta G$ (kJ/mol)',fontsize=f_size)
ax1.set_xlabel(r'$r_{cutoff}$ (nm)',fontsize=f_size)
ax.set_xlabel(r'$r_{cutoff}$ (nm)',fontsize=f_size)

ax.tick_params(axis='both', which='major', labelsize=t_size)
ax1.tick_params(axis='both', which='major', labelsize=t_size)

ax.legend(fontsize=leg_size)
ax1.legend(fontsize=leg_size)

plt.tight_layout()
fig.savefig("dG_vs_rcutoff.svg",format='svg',transparent=True)
plt.show()
