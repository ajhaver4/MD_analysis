import pandas
import sys
import numpy as np
import matplotlib.pyplot as plt

data = pandas.read_csv(sys.argv[1],comment='#',delimiter='\t')
states_data = pandas.read_csv(sys.argv[2],delimiter='\s+')
states={}
for i in range(len(states_data['st'])):

    states[states_data['st'][i]] = 0.0

#Membrane
V = 37.88*37.88*18.03   #nm3
A = 37.88*37.88   #nm2
c0 = 1e3 * (1/1660)    #nm3

# std_c = -8.314*310*np.log(V/c0)/1000    #kJ/mol
std_c=0

print("Standard State Correction :", std_c)
# sys.exit()
print("Timestamp values of dG: ")
print("Time (us) \t\t dG (kJ/mol)")
# times_print = [30,60,90,120]  #For pseudo membrane
# times_print = [90,236]   #For membrane
times_print = [50,130,264]  # For martini2

avg_dG = []
for t in times_print:
    time_val = data['Timestep'].to_numpy()
    indx = np.argmin(np.abs(time_val - t*1e6))
    
    print(f"{data['Timestep'][indx]} \t\t {data['TS1'][indx]}")
    avg_dG.append(data['TS1'][indx])

print("Averaged dG: ", np.mean(avg_dG))
print("SEM: ", np.std(avg_dG)/np.sqrt(len(avg_dG)))    
fig,ax=plt.subplots()

color_st = {'CRY':'black','NS1':'magenta','NS2':'olive','NS3':'steelblue','NS4':'orange','TS2':'royalblue'}
label_st = {'CRY':'Specific','NS1':'Non-Specific(1)','NS2':'Non-Specific(2)','NS3':'Non-Specific(3)','NS4':'Non-Specific(4)','TS2':r'$\Delta G_{2D} $'}

for st,v in states.items():
    if ('unb' not in st) and ('TS1' not in st):
    # if ('TS2' in st) or ('CRY' in st):
    #if ('TS2' in st):
        # print(st,data[st])
        new_data = data[st] > -500
        ax.plot(data['Timestep']/1e6,data[st],linewidth=2.5,alpha=0.8,label=label_st[st],color=color_st[st])
        # ax.plot(data['Timestep'][new_data]/1e6,data[st][new_data],linewidth=4.0,alpha=0.8,label=label_st[st],color=color_st[st])
        print("State: ", st)
        free_energy = data[st].to_numpy()[-1]
        free_energy_std = data[st].to_numpy()[-1] + std_c

        print("Free Energy : ", free_energy)
        print("Standard Free Energy : ", free_energy_std)

        kd = 1e9*np.exp((free_energy_std)*1000/(8.314*310))
        #kd = 1e9*np.exp((data[st].to_numpy()[-1])*1000/(8.314*310))
        # kd=0
        k2d = (1/A)*np.exp((free_energy_std)*1000/(8.314*310))
        print("State: ", st)
        print("Binding Affinity : %.2f nM \t\t dG: %.2f" %(kd,free_energy))
        print(f"2D binding Affinity: %.3E 1/nm\N{SUPERSCRIPT TWO}" %(k2d))


# ax.hlines(np.mean(avg_dG),0,data['Timestep'].to_numpy()[-1]/1e6,linestyle='--',color='black',linewidth=1.5)
ax.legend(loc='best',fontsize=20)
ax.set_ylabel(r'$\Delta G$ kJ/mol',fontsize=30)
ax.set_xlabel(r"Time in $\mu s$",fontsize=30)
ax.tick_params(labelsize=20)
plt.show()
