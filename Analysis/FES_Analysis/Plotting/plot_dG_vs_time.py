import pandas
import sys
import numpy as np
import matplotlib.pyplot as plt

data = pandas.read_csv(sys.argv[1],comment='#',delimiter='\t')
states_data = pandas.read_csv(sys.argv[2],delimiter='\s+')
states={}
print(data)
for i in range(len(states_data['st'])):

    states[states_data['st'][i]] = 0.0

#Membrane

V = 37.88*37.88*37.83   #nm3
c0 = 1e3 * (1/1660)    #nm3

std_c = -8.314*310*np.log(V/c0)/1000    #kJ/mol
# std_c=0

print("Standard State Correction :", std_c)
print("Timestamp values of dG: ")
print("Time (us) \t\t dG (kJ/mol)")
times_print = [30,73,106,165,208]
avg_dG = []
for t in times_print:
    time_val = data['Timestep'].to_numpy()
    indx = np.argmin(np.abs(time_val - t*1e6))
    
    print(f"{data['Timestep'][indx]} \t\t {data['TS1'][indx]}")
    avg_dG.append(data['TS1'][indx])

print("Averaged dG: ", np.mean(avg_dG))
print("SEM: ", np.std(avg_dG)/np.sqrt(len(avg_dG)))    
# sys.exit()
fig_size=(6,4)
f_size=20
t_size=15
leg_size = 14
lw=2.0
fig,ax = plt.subplots(figsize=fig_size)


color_st = {'CRY':'black','NS1':'magenta','NS2':'olive','NS3':'darkorange','NS4':'orange','TS2':'royalblue'}
label_st = {'CRY':'Specific','NS1':'NS-1','NS2':'NS-3','NS3':'NS-2','NS4':'Non-Specific(4)','TS2':r'$\Delta G_{3D} $'}

min_dGs = []
for st,v in states.items():
    # if ('unb' not in st) and ('TS2' not in st):
    if ('TS2' in st):
        # print(st,data[st])
        # new_data = data[st] > -500
        ax.plot(data['Timestep']/1e6,data[st],linewidth=lw,alpha=0.8,label=label_st[st],color=color_st[st])
        # ax.plot(data['Timestep'][new_data]/1e6,data[st][new_data],linewidth=4.0,alpha=0.8,label=label_st[st],color=color_st[st])
        print("State: ", st)
        free_energy = data[st].to_numpy()[-1]
        free_energy_std = data[st].to_numpy()[-1] + std_c

        print("Free Energy : ", free_energy)
        print("Standard Free Energy : ", free_energy_std)

        kd = 1e9*np.exp((free_energy_std)*1000/(8.314*310))
        #kd = 1e9*np.exp((data[st].to_numpy()[-1])*1000/(8.314*310))
        # kd=0
        k3d = (1/V)*np.exp((free_energy)*1000/(8.314*310))
        print("State: ", st)
        print("Binding Affinity (Using std conc.): %.2f nM \t\t dG: %.2f" %(kd,free_energy))
        print(f"3D binding Affinity: %.3E 1/nm\N{SUPERSCRIPT THREE}" %(k3d))

        min_dGs.append(np.min(data[st].to_numpy()))

# ax.hlines(np.mean(avg_dG),0,data['Timestep'].to_numpy()[-1]/1e6,linestyle='--',color='black',linewidth=1.5)

total_time = data['Timestep'].to_numpy()[-1]/1e6
half_time = total_time/2
block_times = np.arange(half_time,total_time,20)

min_dG = np.min(min_dGs)
start_time=half_time
ax.axvspan(0,half_time,facecolor='gray',alpha=0.3)
# for i in range(1,len(block_times)):
#     end_time = block_times[i]
#     if i == len(block_times)-1:
#         end_time = total_time
#     ax.axvspan(start_time,end_time,facecolor='crimson',alpha=0.3)
#     ax.vlines(start_time,min_dG-1,0,linestyle='--',color='black',linewidth=1.2)
#     start_time = block_times[i]

ax.axvspan(half_time,half_time+20,facecolor='crimson',alpha=0.3)


ax.legend(loc='lower right',fontsize=leg_size)
ax.set_ylabel(r'$\Delta G$ kJ/mol',fontsize=f_size)
ax.set_xlabel(r"Time in $\mu s$",fontsize=f_size)
ax.tick_params(labelsize=t_size)
ax.set_ylim(bottom=min_dG-1,top=0)
fig.tight_layout()
plt.savefig("dG_vs_time.svg",dpi=600,format="svg",transparent=True)
plt.show()



