import pandas
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

import numpy as np
from scipy.stats import ks_2samp

data3D = pandas.read_csv(sys.argv[1],comment='#',delimiter='\t')
data2D = pandas.read_csv(sys.argv[2],comment='#',delimiter='\t')
# dataPs = pandas.read_csv(sys.argv[3],comment='#',delimiter='\t')


#3D, 2D and 2D-Ps
avg_start_times = [60,100]
# avg_dG_values = [-23.341,-95.972,-49.434]

times3D = data3D['Timestep'].to_numpy()
times2D = data2D['Timestep'].to_numpy()
# timesPs = data3D['Timestep'].to_numpy()

indx3D = np.argmin(np.abs(times3D - avg_start_times[0]*1e6))
indx2D = np.argmin(np.abs(times2D - avg_start_times[1]*1e6))
# indxPs = np.argmin(np.abs(timesPs - avg_start_times[2]*1e6))

avg_data3D = data3D[indx3D:]
avg_data2D = data2D[indx2D:]
# avg_dataPs = dataPs[indxPs:]

print(avg_data2D['Timestep'])


ts_3D = times3D[1]-times3D[0]
ts_2D = times2D[1]-times2D[0]
# ts_PS = timesPs[1]-timesPs[0]

std_timestamps = np.arange(1,1.1,0.5)   #us
std_timestamps3D_ctrl = np.arange(2,50,0.7)
std_timestamps2D = np.arange(1,1.1,4)

time_int3D = np.array([int(i*1e6/ts_3D) for i in std_timestamps])
time_int3D_ctrl = np.array([int(i*1e6/ts_3D) for i in std_timestamps3D_ctrl])
time_int2D = np.array([int(i*1e6/ts_2D) for i in std_timestamps2D])
# time_intPs = np.array([int(i*1e6/ts_PS) for i in std_timestamps])

#Get the dG's for each time interval
dG_blocks3D = {std_timestamps[i] : avg_data3D['TS2'][::time_int3D[i]] for i in range(len(std_timestamps))}
dG_blocks3D_ctrl = {std_timestamps3D_ctrl[i] : avg_data3D['TS2'][::time_int3D_ctrl[i]] for i in range(len(std_timestamps3D_ctrl))}
dG_blocks2D = {std_timestamps2D[i] : avg_data2D['TS2'][::time_int2D[i]] for i in range(len(std_timestamps2D))}
# dG_blocksPs = {std_timestamps[i] : avg_dataPs['TS2'][::time_intPs[i]] for i in range(len(std_timestamps))}

std_values3D = {std_timestamps[i] : np.std(dG_blocks3D[std_timestamps[i]],ddof=1) for i in range(len(std_timestamps))}
std_values3D_ctrl = {std_timestamps3D_ctrl[i] : np.std(dG_blocks3D_ctrl[std_timestamps3D_ctrl[i]],ddof=1) for i in range(len(std_timestamps3D_ctrl))}
std_values2D = {std_timestamps2D[i] : np.std(dG_blocks2D[std_timestamps2D[i]],ddof=1) for i in range(len(std_timestamps2D))}
# std_valuesPs = {std_timestamps[i] : np.std(dG_blocksPs[std_timestamps[i]],ddof=1) for i in range(len(std_timestamps))}

avg_values3D = {std_timestamps[i] : np.mean(dG_blocks3D[std_timestamps[i]]) for i in range(len(std_timestamps))}
avg_values3D_ctrl = {std_timestamps3D_ctrl[i] : np.mean(dG_blocks3D_ctrl[std_timestamps3D_ctrl[i]]) for i in range(len(std_timestamps3D_ctrl))}
avg_values2D = {std_timestamps2D[i] : np.mean(dG_blocks2D[std_timestamps2D[i]]) for i in range(len(std_timestamps2D))}
# avg_valuesPs = {std_timestamps[i] : np.mean(dG_blocksPs[std_timestamps[i]]) for i in range(len(std_timestamps))}

cf = 0.95    # % CI
norm_error = norm.ppf((cf+1)/2)
error_bars3D = { i : norm_error*std_values3D[i]/np.sqrt(len(dG_blocks3D[i])) for i in std_timestamps}
error_bars3D_ctrl = { i : norm_error*std_values3D_ctrl[i]/np.sqrt(len(dG_blocks3D_ctrl[i])) for i in std_timestamps3D_ctrl}
error_bars2D = { i : norm_error*std_values2D[i]/np.sqrt(len(dG_blocks2D[i])) for i in std_timestamps2D}
# error_barsPs = { i : norm_error*std_valuesPs[i]/np.sqrt(len(dG_blocksPs[i])) for i in std_timestamps}


hist_values3D = [val for dG_values in dG_blocks3D.values() for val in dG_values]
hist_values3D_ctrl = [val for dG_values in dG_blocks3D_ctrl.values() for val in dG_values]
hist_values2D = [val for dG_values in dG_blocks2D.values() for val in dG_values]


#Shift 3D values by the required ddG so that the resulatant distribution can be compared to the actual 2D distribution for hypothesis testing
hyp_estimate2D = [val - 23.52 for val in hist_values3D] 
print(len(hist_values3D))
print(len(hist_values2D))


def get_pvalues(estimated_values,actual_values,ax=None,label=None,color=None):
    # ax.hist(hist_values3D,bins=10,label='3D',alpha=0.5,color='navy')
    ax.hist(estimated_values,bins=40,label=r'$\Delta G_{rigid}$',alpha=0.5,color='gray',density=True)
    ax.hist(actual_values,bins=20,label=label,alpha=0.5,color=color,density=True)

    ax.legend(fontsize=21,loc='upper left',frameon=False)
    ax.set_xlabel(r'$\Delta G (kJ/mol)$',fontsize=26)
    ax.set_ylabel(r'P $(\Delta G)$',fontsize=26)
    ax.tick_params(labelsize=24)


    # Perform KS test
    ks_stat, p_value = ks_2samp(estimated_values, actual_values)

    # ax.text(0.05, 0.5, f"p-value: {p_value:.3e}", transform=ax.transAxes,fontsize=16)

    return ks_stat, p_value




fig,ax = plt.subplots(figsize=(6,4))

ks_stat,p_value = get_pvalues(hyp_estimate2D,hist_values2D,ax=ax,label=r'$\Delta G_{2D}$',color='firebrick')
fig.tight_layout()
fig.savefig('hvalues_test_2DMart2.png',dpi=600,transparent=True,format='png')
fig3,ax3 = plt.subplots(figsize=(6,4))

ks_stat3D,p_value3D = get_pvalues(hist_values3D,hist_values3D_ctrl,ax=ax3,label='3D control',color='navy')
fig3.tight_layout()

plt.show()

with open('dG_FES_3D', 'w') as fl:
    fl.write('#dG_rigid\n')
    for val in hist_values3D:
        fl.write(f'{val}\n')
        
# with open('dG_Rigid_2D_M2', 'w') as fl:
#     fl.write('#dG_rigid\n')
#     for val in hyp_estimate2D:
#         fl.write(f'{val}\n')

# with open('dG_FES_2D_M2', 'w') as fl:
#     fl.write('#dG_2D\n')
#     for val in hist_values2D:
#         fl.write(f'{val}\n')
print("-------------------------------------------------------------")

#Interpret the result
alpha = 0.05

print("For the 2D and 3D distributions: ")
print(f"KS Statistic: {ks_stat}")
print(f"P-value: {p_value}")

if p_value > alpha:
    print("The null hypothesis cannot be rejected. The distributions are likely from the same distribution.")
else:
    print("The null hypothesis can be rejected. The measurements are likely from different distributions which indicates no significant overlap")

print("-------------------------------------------------------------")
print("For the 3D and 3D control distributions: ")
print(f"KS Statistic 3D: {ks_stat3D}")
print(f"P-value 3D: {p_value3D}")
if p_value3D > alpha:
    print("The null hypothesis cannot be rejected. The distributions are likely from the same distribution.")
else:
    print("The null hypothesis can be rejected. The measurements are likely from different distributions which indicates no significant overlap")



"""
Checking with only Averaged values
"""
#Shift 3D values by the required ddG so that the resulatant distribution can be compared to the actual 2D distribution for hypothesis testing
hyp_estimate2D = [val - 23.52 for val in avg_values3D.values()] 

fig,ax = plt.subplots(figsize=(10,6))
ks_stat,p_value = get_pvalues(hyp_estimate2D,list(avg_values2D.values()),ax=ax,label=r'$\Delta G_{2D}$',color='firebrick')

fig3,ax3 = plt.subplots(figsize=(10,6))
ks_stat3D,p_value3D = get_pvalues(list(avg_values3D.values()),list(avg_values3D_ctrl.values()),ax=ax3,label='3D control',color='navy')

plt.show()

print("-------------------------------------------------------------")
#Interpret the result
alpha = 0.05
print("For the 2D and 3D distributions: ")
print(f"KS Statistic: {ks_stat}")
print(f"P-value: {p_value}")
if p_value > alpha:
    print("The null hypothesis cannot be rejected. The distributions are likely from the same distribution.")
else:
    print("The null hypothesis can be rejected. The measurements are likely from different distributions which indicates no significant overlap")

print("-------------------------------------------------------------")
print("For the 3D and 3D control distributions: ")
print(f"KS Statistic 3D: {ks_stat3D}")
print(f"P-value 3D: {p_value3D}")
if p_value3D > alpha:
    print("The null hypothesis cannot be rejected. The distributions are likely from the same distribution.")
else:
    print("The null hypothesis can be rejected. The measurements are likely from different distributions which indicates no significant overlap")


