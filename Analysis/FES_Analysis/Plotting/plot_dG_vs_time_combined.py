import pandas
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

data3D = pandas.read_csv(sys.argv[1],comment='#',delimiter='\t')
data2D = pandas.read_csv(sys.argv[2],comment='#',delimiter='\t')
dataPs = pandas.read_csv(sys.argv[3],comment='#',delimiter='\t')


#3D, 2D and 2D-Ps
avg_start_times = [104,145,85]
# avg_dG_values = [-23.341,-95.972,-49.434]

times3D = data3D['Timestep'].to_numpy()
times2D = data3D['Timestep'].to_numpy()
timesPs = data3D['Timestep'].to_numpy()

indx3D = np.argmin(np.abs(times3D - avg_start_times[0]*1e6))
indx2D = np.argmin(np.abs(times2D - avg_start_times[1]*1e6))
indxPs = np.argmin(np.abs(timesPs - avg_start_times[2]*1e6))

avg_data3D = data3D[indx3D:]
avg_data2D = data2D[indx2D:]
avg_dataPs = dataPs[indxPs:]



ts_3D = times3D[1]-times3D[0]
ts_2D = times2D[1]-times2D[0]
ts_PS = timesPs[1]-timesPs[0]

std_timestamps = [0.1,0.5,1,1.5,2,5,10,15,20,25,30,35,40,45,50]   #us

time_int3D = np.array([int(i*1e6/ts_3D) for i in std_timestamps])
time_int2D = np.array([int(i*1e6/ts_2D) for i in std_timestamps])
time_intPs = np.array([int(i*1e6/ts_PS) for i in std_timestamps])

#Get the dG's for each time interval
dG_blocks3D = {std_timestamps[i] : avg_data3D['TS2'][::time_int3D[i]] for i in range(len(std_timestamps))}
dG_blocks2D = {std_timestamps[i] : avg_data2D['TS2'][::time_int2D[i]] for i in range(len(std_timestamps))}
dG_blocksPs = {std_timestamps[i] : avg_dataPs['TS2'][::time_intPs[i]] for i in range(len(std_timestamps))}

std_values3D = {std_timestamps[i] : np.std(dG_blocks3D[std_timestamps[i]],ddof=1) for i in range(len(std_timestamps))}
std_values2D = {std_timestamps[i] : np.std(dG_blocks2D[std_timestamps[i]],ddof=1) for i in range(len(std_timestamps))}
std_valuesPs = {std_timestamps[i] : np.std(dG_blocksPs[std_timestamps[i]],ddof=1) for i in range(len(std_timestamps))}

avg_values3D = {std_timestamps[i] : np.mean(dG_blocks3D[std_timestamps[i]]) for i in range(len(std_timestamps))}
avg_values2D = {std_timestamps[i] : np.mean(dG_blocks2D[std_timestamps[i]]) for i in range(len(std_timestamps))}
avg_valuesPs = {std_timestamps[i] : np.mean(dG_blocksPs[std_timestamps[i]]) for i in range(len(std_timestamps))}

cf = 0.95    # % CI
norm_error = norm.ppf((cf+1)/2)
error_bars3D = { i : norm_error*std_values3D[i]/np.sqrt(len(dG_blocks3D[i])) for i in std_timestamps}
error_bars2D = { i : norm_error*std_values2D[i]/np.sqrt(len(dG_blocks2D[i])) for i in std_timestamps}
error_barsPs = { i : norm_error*std_valuesPs[i]/np.sqrt(len(dG_blocksPs[i])) for i in std_timestamps}



#For the plot line select which block size we wnat to use to plot the corresponding avg dG 
avg_dGs = [avg_values3D[20],avg_values2D[40],avg_valuesPs[30]]



std_data3D = np.std(avg_data3D['TS2'])
std_data2D = np.std(avg_data2D['TS2'])
std_dataPS = np.std(avg_dataPs['TS2'])

print("Standard Deviation 3D: ", std_data3D)
print("Standard Deviation 2D: ", std_data2D)
print("Standard Deviation 2D-Ps: ", std_dataPS)

print("--------------------------------------------------------")
print("Average Values 3D: ", avg_values3D)
print("Average Values 2D: ", avg_values2D)
print("Average Values 2D-Ps: ", avg_valuesPs)

print("----------------------------------------------------------")
print("Standard Deviation for a range of different time intervals: ")
print("Standard Deviation 3D: ", std_values3D)
print("Standard Deviation 2D: ", std_values2D)
print("Standard Deviation 2D-Ps: ", std_valuesPs)

print("-------------------------------------------------------------")
print("Error Bars for a range of different time intervals: ")
print("Error Bars 3D: ", error_bars3D)
print("Error Bars 2D: ", error_bars2D)
print("Error Bars 2D-Ps: ", error_barsPs)

fig_size=(12,4)
f_size=24
t_size=18
leg_size = 18
lw=2.0
fig,[ax,ax1,ax2] = plt.subplots(1,3,figsize=fig_size,sharey=True)

ax.plot(data3D['Timestep']/1e6,data3D['TS2'],linewidth=lw,alpha=0.8,label=r'$\Delta G_{3D} $',color='black')
ax1.plot(data2D['Timestep']/1e6,data2D['TS2'],linewidth=lw,alpha=0.8,label=r'$\Delta G_{2D} $',color='black')
ax2.plot(dataPs['Timestep']/1e6,dataPs['TS2'],linewidth=lw,alpha=0.8,label=r'$\Delta G_{2D-Ps} $',color='black')

ax.hlines(avg_dGs[0], 0, data3D['Timestep'].to_numpy()[-1]/1e6, colors='black', linestyles='dashed')
ax1.hlines(avg_dGs[1], 0, data2D['Timestep'].to_numpy()[-1]/1e6, colors='black', linestyles='dashed')
ax2.hlines(avg_dGs[2], 0, dataPs['Timestep'].to_numpy()[-1]/1e6, colors='black', linestyles='dashed')



ax.set_ylabel(r'$\Delta G$ (kJ/mol)',fontsize=f_size)
ax.set_xlabel('Time (us)',fontsize=f_size)
ax1.set_xlabel('Time (us)',fontsize=f_size)
ax2.set_xlabel('Time (us)',fontsize=f_size)

ax.tick_params(axis='both', which='major', labelsize=t_size)
ax1.tick_params(axis='both', which='major', labelsize=t_size)
ax2.tick_params(axis='both', which='major', labelsize=t_size)

ax.legend(fontsize=leg_size)
ax1.legend(fontsize=leg_size)
ax2.legend(fontsize=leg_size)

plt.tight_layout()
fig.savefig("dG_vs_time.svg",format='svg',transparent=True)
plt.show()


fig,[ax,ax1,ax2] = plt.subplots(1,3,figsize=fig_size,sharey=True)

ax.plot(std_timestamps,avg_values3D.values(),linewidth=lw,alpha=0.8,label=r'$\Delta G_{3D} $',color='black')
ax.errorbar(std_timestamps,avg_values3D.values(),yerr=error_bars3D.values(),fmt='o',color='black',capsize=5)

ax1.plot(std_timestamps,avg_values2D.values(),linewidth=lw,alpha=0.8,label=r'$\Delta G_{2D} $',color='black')
ax1.errorbar(std_timestamps,avg_values2D.values(),yerr=error_bars2D.values(),fmt='o',color='black',capsize=5)

ax2.plot(std_timestamps,avg_valuesPs.values(),linewidth=lw,alpha=0.8,label=r'$\Delta G_{2D-Ps} $',color='black')
ax2.errorbar(std_timestamps,avg_valuesPs.values(),yerr=error_barsPs.values(),fmt='o',color='black',capsize=5)

ax.set_ylabel(r'$\Delta G$ (kJ/mol)',fontsize=f_size)
ax.set_xlabel('Block Size (us)',fontsize=f_size)
ax1.set_xlabel('Block Size (us)',fontsize=f_size)
ax2.set_xlabel('Block Size (us)',fontsize=f_size)

ax.tick_params(axis='both', which='major', labelsize=t_size)
ax1.tick_params(axis='both', which='major', labelsize=t_size)
ax2.tick_params(axis='both', which='major', labelsize=t_size)

ax.set_xticks([0,10,20,30,40,50])
ax1.set_xticks([0,10,20,30,40,50])
ax2.set_xticks([0,10,20,30,40,50])

ax.legend(fontsize=leg_size)
ax1.legend(fontsize=leg_size)
ax2.legend(fontsize=leg_size)


plt.tight_layout()
fig.savefig("dG_error.svg",format='svg',transparent=True)
plt.show()




with open("dG_blockData3D.dat",'w') as fl:
    fl.write("#Block Time Interval (us) \t Average dG (kJ/mol) \t STD \t Error")
    fl.write("\n")
    for i in std_timestamps:
        fl.write("%2.2f \t %2.2f \t %2.2f \t %2.2f" %(i,avg_values3D[i],std_values3D[i],error_bars3D[i]))
        fl.write("\n")

    fl.write("%2.2f \t %2.2f \t %2.2f \t %2.2f" %(avg_data3D['Timestep'].values[-1],avg_data3D['TS2'].values[-1],0.0,0.0))
    fl.write("\n")


with open("dG_blockData2D.dat",'w') as fl:
    fl.write("#Block Time Interval (us) \t Average dG (kJ/mol) \t STD \t Error")
    fl.write("\n")
    for i in std_timestamps:
        fl.write("%2.2f \t %2.2f \t %2.2f \t %2.2f" %(i,avg_values2D[i],std_values2D[i],error_bars2D[i]))
        fl.write("\n")

    fl.write("%2.2f \t %2.2f \t %2.2f \t %2.2f" %(avg_data2D['Timestep'].values[-1],avg_data2D['TS2'].values[-1],0.0,0.0))
    fl.write("\n")


with open("dG_blockDataPs.dat",'w') as fl:
    fl.write("#Block Time Interval (us) \t Average dG (kJ/mol) \t STD \t Error")
    fl.write("\n")
    for i in std_timestamps:
        fl.write("%2.2f \t %2.2f \t %2.2f \t %2.2f" %(i,avg_valuesPs[i],std_valuesPs[i],error_barsPs[i]))
        fl.write("\n")

    fl.write("%2.2f \t %2.2f \t %2.2f \t %2.2f" %(avg_dataPs['Timestep'].values[-1],avg_dataPs['TS2'].values[-1],0.0,0.0))
    fl.write("\n")


