import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from matplotlib.markers import MarkerStyle as markerstyle
import math
from scipy.interpolate import griddata
import pandas
from scipy import constants
import matplotlib.patches as patches
import matplotlib.cm as cm
from Clust_times import time_dict

"""
Inputs - Final_data fes_file energy_file
"""

"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-------------------------------- INITIALIZATION --------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
"""
#Reading in the structure file and trajectory
dist_data = pandas.read_csv(sys.argv[1],comment='#',sep='\t')    #Master distances file
fes_file = sys.argv[2]
ener_data = pandas.read_csv(sys.argv[3],comment='#',sep='\t')   #Reweighted Energy file

"""
------- Defining Constants ------------
"""
k = constants.value('Boltzmann constant')
Ava_no = constants.value('Avogadro constant')
temp = 310
kbt_ = 1000/(k*temp*Ava_no)

"""
------------- Defining Useful methods ----------------
"""
def get_indx_from_df(df, col_name, val):
   mod_val = df[col_name]-val
   mod_val = mod_val.abs()
   return(mod_val.idxmin())
#    return (df.index[df[col_name] == val].tolist())


def get_free_energy(x_cood,y_cood,X_grids,Y_grids):
    diff = Y_grids-y_cood
    bool = diff<=0
    bool_arr = diff[bool]
    max = np.amax(bool_arr)
    idx_y = list(diff).index(max)

    diff = X_grids-x_cood
    bool = diff<=0
    bool_arr = diff[bool]
    max = np.amax(bool_arr)
    idx_x = list(diff).index(max)

    return((idx_x,idx_y))



"""
--------------------------------------------------------------------------------
--------------------------------- READING FES DATA -----------------------------
--------------------------------------------------------------------------------
"""

print("Reading FES data...")
data_fes = [line.split() for line in open(fes_file, "r")]
data2 = [x for x in data_fes if not x == []]  # strip out headers
# file_name = str(sys.argv[5])
d1, d2, free, dd1, dd2 = [], [], [], [], []
for elem in data2[9:]:
    d1.append(float(elem[0]))
    d2.append(float(elem[1]))
    free.append(float(elem[2]))
#    dd1.append(float(elem[3]))
#    dd2.append(float(elem[4]))

X = np.linspace(min(d1), max(d1), 1318)
Y = np.linspace(min(d2), max(d2), 1322)

#Normalize unbound state to zero
free = np.array(free)
# free = free+110.6

#Zero value
d1_arr = np.array(d1)
d2_arr = np.array(d2)
mask1 = (d1_arr >= 25.0 ) & (d1_arr < 26.0 )
mask2 = (d2_arr >= 25.0 ) & (d2_arr < 26.0 )
mask3 = mask1 * mask2

correction = np.mean(free[mask3])
free = free - correction

#Shift max value to a constant. To shifht the color map
max_val = 50.0
mask4 = (free>=max_val)
free[mask4]=max_val

print("Creating data grid. This may take a while...")
D1, D2 = np.meshgrid(X, Y)
ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)


"""
--------------------------------------------------------------------------------
--------------------------------- READING CLUSTER DATA -------------------------
--------------------------------------------------------------------------------
"""

cluster_states = list(set(time_dict.values()))
cluster_states = cluster_states +['Bound','Unb12','Unb13','Unb14']

cluster_thermo = {st:{'potenergy':[],'rbias':[],'free':[],'volume':[],'d1':[],'d2':[]} for st in cluster_states}
cluster_thermo_df = pandas.DataFrame(columns=['Cluster','Bonded','Col-Pro','LJ-Pro','Col-ProMem','LJ-ProMem','Col-ProW','LJ-ProW','Col-W','LJ-W','Col-Mem','LJ-Mem','PotEnergy','FreeEnergy','Volume','rbias','d1','d2'])

print("Reading REWEIGHTED data...")

time_errors = []
# for time,states in time_dict.items():
#     try:
#         time_indx = get_indx_from_df(ener_data,'Timestep',time)
#         # time_indx = time_indx[0]
#     except Exception as e:
#         print("Error: ",e)
#         print("Time not found in energy file: ",time,time_indx)
#         time_errors.append(time)
#         continue

#     cluster_thermo[states]['potenergy'].append(ener_data['Pot'][time_indx])
#     cluster_thermo[states]['rbias'].append(ener_data['rbias'][time_indx])

    
#     bonded_term = ener_data['Bond'][time_indx]+ener_data['HP'][time_indx]+ener_data['G96'][time_indx]+ener_data['PDih'][time_indx]+ener_data['ImPDih'][time_indx]
#     pro_col = ener_data['Col-Pro-Pro'][time_indx]
#     proMem_col = ener_data['Col-Pro-Mem'][time_indx]
#     proW_col = ener_data['Col-Pro-W'][time_indx]
#     pro_LJ = ener_data['LJ-Pro-Pro'][time_indx]
#     proMem_LJ = ener_data['LJ-Pro-Mem'][time_indx]
#     proW_LJ = ener_data['LJ-Pro-W'][time_indx]
#     sol_col = ener_data['Col-W-W'][time_indx]
#     sol_LJ = ener_data['LJ-W-W'][time_indx]
#     mem_col = ener_data['Col-Mem-Mem'][time_indx]+ener_data['Col-Mem-W'][time_indx]
#     mem_LJ = ener_data['LJ-Mem-Mem'][time_indx]+ener_data['LJ-Mem-W'][time_indx]
    

#     g = get_free_energy(ener_data['d1'][time_indx],ener_data['d2'][time_indx],X,Y)
#     free_energy = ENER[g[1]][g[0]]
#     cluster_thermo[states]['free'].append(free_energy)

#     cluster_thermo[states]['volume'].append(ener_data['Volume'][time_indx])
#     cluster_thermo[states]['d1'].append(ener_data['d1'][time_indx])
#     cluster_thermo[states]['d2'].append(ener_data['d2'][time_indx])

#     new_data = [states,bonded_term,pro_col,pro_LJ,proMem_col,proMem_LJ,proW_col,proW_LJ,sol_col,sol_LJ,mem_col,mem_LJ,ener_data['Pot'][time_indx],free_energy,ener_data['Volume'][time_indx],
#                 ener_data['rbias'][time_indx],ener_data['d1'][time_indx],ener_data['d2'][time_indx]]
#     cluster_thermo_df = pandas.concat([pandas.DataFrame([new_data],columns=cluster_thermo_df.columns),cluster_thermo_df],ignore_index=True)



final_cluster_stats = {st:{'potenergy':[],'free':[],'volume':[]} for st in cluster_states}


#Calculating for unbound state
unb1 = 12.0
unb2 = 13.0
unb3 = 14.0



#Filtering of unbound states needs to be done on a frame by frame basis since the timestamps of 
# MDAnalysis and gromacs energy output do not match exactly (possibly due ot round off error). 
for i in range(len(dist_data['Timestep'])):
    unb_labels=[]
    if dist_data['d1'][i] <= 19.0 and dist_data['d2'][i]<=19.0:
        if dist_data['d1'][i] >= unb1 and dist_data['d2'][i]>=unb1:
            unb_labels.append('Unb12')
            try:
                time = dist_data['Timestep'][i]
                time_indx = get_indx_from_df(ener_data,'Timestep',time)
                time_indx = time_indx
            except Exception as e:
                print("Error: ",e)
                print("Time not found in energy file: ",time,time_indx)
                time_errors.append(time)
                continue
            
        if dist_data['d1'][i]>=unb2 and dist_data['d2'][i]>=unb2:
            unb_labels.append('Unb13')
        if dist_data['d1'][i]>=unb3 and dist_data['d2'][i]>=unb3:
            unb_labels.append('Unb14')
        
        for ulabel in unb_labels: 
            cluster_thermo[ulabel]['potenergy'].append(ener_data['Pot'][time_indx])
            cluster_thermo[ulabel]['rbias'].append(ener_data['rbias'][time_indx])

            g = get_free_energy(ener_data['d1'][time_indx],ener_data['d2'][time_indx],X,Y)
            cluster_thermo[ulabel]['free'].append(ENER[g[1]][g[0]])
            cluster_thermo[ulabel]['volume'].append(ener_data['Volume'][time_indx])
            cluster_thermo[ulabel]['d1'].append(ener_data['d1'][time_indx])
            cluster_thermo[ulabel]['d2'].append(ener_data['d2'][time_indx])

            #Extracting all the Energy terms for the unbound state
            bonded_term = ener_data['Bond'][time_indx]+ener_data['HP'][time_indx]+ener_data['G96'][time_indx]+ener_data['PDih'][time_indx]+ener_data['ImPDih'][time_indx]
            pro_col = ener_data['Col-Pro-Pro'][time_indx]
            proMem_col = ener_data['Col-Pro-Mem'][time_indx]
            proW_col = ener_data['Col-Pro-W'][time_indx]
            pro_LJ = ener_data['LJ-Pro-Pro'][time_indx]
            proMem_LJ = ener_data['LJ-Pro-Mem'][time_indx]
            proW_LJ = ener_data['LJ-Pro-W'][time_indx]
            sol_col = ener_data['Col-W-W'][time_indx]
            sol_LJ = ener_data['LJ-W-W'][time_indx]
            mem_col = ener_data['Col-Mem-Mem'][time_indx]+ener_data['Col-Mem-W'][time_indx]
            mem_LJ = ener_data['LJ-Mem-Mem'][time_indx]+ener_data['LJ-Mem-W'][time_indx]

            new_data = [ulabel,bonded_term,pro_col,pro_LJ,proMem_col,proMem_LJ,proW_col,proW_LJ,sol_col,sol_LJ,mem_col,mem_LJ,ener_data['Pot'][time_indx],free_energy,ener_data['Volume'][time_indx],
                ener_data['rbias'][time_indx],ener_data['d1'][time_indx],ener_data['d2'][time_indx]]
            cluster_thermo_df = pandas.concat([pandas.DataFrame([new_data],columns=cluster_thermo_df.columns),cluster_thermo_df],ignore_index=True)



        
        if dist_data['d1'][i] < unb1 and dist_data['d2'][i]<unb1:
            try:
                time = dist_data['Timestep'][i]
                time_indx = get_indx_from_df(ener_data,'Timestep',time)
                time_indx = time_indx
            except Exception as e:
                print("Error: ",e)
                print("Time not found in energy file: ",time,time_indx)
                time_errors.append(time)
                continue
            bdlabel='Bound'
            cluster_thermo[bdlabel]['potenergy'].append(ener_data['Pot'][time_indx])
            cluster_thermo[bdlabel]['rbias'].append(ener_data['rbias'][time_indx])

            g = get_free_energy(ener_data['d1'][time_indx],ener_data['d2'][time_indx],X,Y)
            cluster_thermo[bdlabel]['free'].append(ENER[g[1]][g[0]])
            cluster_thermo[bdlabel]['volume'].append(ener_data['Volume'][time_indx])
            cluster_thermo[bdlabel]['d1'].append(ener_data['d1'][time_indx])
            cluster_thermo[bdlabel]['d2'].append(ener_data['d2'][time_indx])

            #Extracting all the Energy terms for the Bound state
            bonded_term = ener_data['Bond'][time_indx]+ener_data['HP'][time_indx]+ener_data['G96'][time_indx]+ener_data['PDih'][time_indx]+ener_data['ImPDih'][time_indx]
            pro_col = ener_data['Col-Pro-Pro'][time_indx]
            proMem_col = ener_data['Col-Pro-Mem'][time_indx]
            proW_col = ener_data['Col-Pro-W'][time_indx]
            pro_LJ = ener_data['LJ-Pro-Pro'][time_indx]
            proMem_LJ = ener_data['LJ-Pro-Mem'][time_indx]
            proW_LJ = ener_data['LJ-Pro-W'][time_indx]
            sol_col = ener_data['Col-W-W'][time_indx]
            sol_LJ = ener_data['LJ-W-W'][time_indx]
            mem_col = ener_data['Col-Mem-Mem'][time_indx]+ener_data['Col-Mem-W'][time_indx]
            mem_LJ = ener_data['LJ-Mem-Mem'][time_indx]+ener_data['LJ-Mem-W'][time_indx]

            new_data = [bdlabel,bonded_term,pro_col,pro_LJ,proMem_col,proMem_LJ,proW_col,proW_LJ,sol_col,sol_LJ,mem_col,mem_LJ,ener_data['Pot'][time_indx],free_energy,ener_data['Volume'][time_indx],
                ener_data['rbias'][time_indx],ener_data['d1'][time_indx],ener_data['d2'][time_indx]]
            cluster_thermo_df = pandas.concat([pandas.DataFrame([new_data],columns=cluster_thermo_df.columns),cluster_thermo_df],ignore_index=True)

    
    
    else:
        continue

# print("Lenght of Unbound frames: ",len(cluster_thermo['Unb12']['potenergy']))
# sys.exit()    
for cluster,data in cluster_thermo.items():
    print("Cluster: ",cluster)
    
    bin_reweights = np.exp(np.array(data['rbias'])*kbt_)
    avg_energy = np.average(np.array(data['potenergy']),weights=bin_reweights,axis=0)
    std_energy = np.std(np.array(data['potenergy']),axis=0)
    nonrew_energy = np.average(np.array(data['potenergy']),axis=0)

    avg_volume = np.average(np.array(data['volume']),axis=0)

    #Integrating FES
    #Integration of FES:
    print("Integrating FES...")
    free_ener = np.array(data['free'])
    free_ener = free_ener[np.nonzero(free_ener)]
    #Defining constants
    k = constants.value('Boltzmann constant')
    Ava_no = constants.value('Avogadro constant')
    temp = 310
    free_ener2 = (free_ener*1000)/(k*temp*Ava_no)
    free_ener3 = np.exp(-1*free_ener2)
    bin1_width = X[1]-X[0]
    bin2_width = Y[1]-Y[0]
    print('Bin widths:', bin1_width,bin2_width)
    free_ener4 = free_ener3*bin1_width*bin2_width
    int_free_ener = np.sum(free_ener4,axis=0)
    final_free = -k*Ava_no*temp*(math.log(int_free_ener))/1000

    final_cluster_stats[cluster]['potenergy'] = avg_energy
    final_cluster_stats[cluster]['free'] = final_free
    final_cluster_stats[cluster]['volume'] = avg_volume

    print("------------------------")
    print("Cluster: ",cluster)
    print("Average Energy: ",avg_energy)
    print("STD Energy: ",std_energy)
    print("Non-rewighted Energy: ",nonrew_energy)
    print("Free Energy: ",final_free)
    print("Volume: ",avg_volume)
    print("------------------------")


    filter_df=cluster_thermo_df[cluster_thermo_df['Cluster']==cluster]
    print("Classified Energy Stats: ")
    print("%-15s%-15s%-15s" %( "Quantity","Mean","Std"))
    columes = ['Bonded','Col-Pro','LJ-Pro','Col-ProMem','LJ-ProMem','Col-ProW','LJ-ProW','Col-W','LJ-W','Col-Mem','LJ-Mem','PotEnergy']
    
    bin_reweights = np.exp(np.array(filter_df['rbias'])*kbt_)
    for i in range(len(columes)):
        print("%-15s" %(columes[i]),end='')
        print("%-13.2f" %(np.average(filter_df[columes[i]],weights=bin_reweights)),end='')
        print("%-13.2f" %(np.std(filter_df[columes[i]])),end='')
        print("\n")



