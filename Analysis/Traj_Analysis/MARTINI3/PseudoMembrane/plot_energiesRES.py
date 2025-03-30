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
from matplotlib.ticker import FormatStrFormatter


ener_data = pandas.read_csv(sys.argv[1],comment='#',sep='\s+',dtype=np.float64)   #Reweighted Energy file
dist_data = pandas.read_csv(sys.argv[2],comment='#',sep='\t',dtype=np.float64)   #Distance file
dist_data['Timestep'] = 200*round(dist_data['Timestep']/200)
#Remove duplicates 
dist_data = dist_data.drop_duplicates(subset=['Timestep'],ignore_index=True)
ener_data = ener_data.drop_duplicates(subset=['Timestep'],ignore_index=True)

# print("Energy Data: ", ener_data['Timestep'])

# print("Final_data: ", dist_data['Timestep'])

# # t_values = {'Timestep': dist_data['Timestep']}.values()
# mask1 = ener_data['Timestep'].isin(dist_data['Timestep'].values)
# print(mask1)
# ener_data = ener_data[~mask1]
# print(ener_data['Timestep'])



master_data = dist_data.join(ener_data.set_index('Timestep'),on='Timestep',how='inner',lsuffix = '_mda',rsuffix = '_plumed')


col_names=['Timestep','Bond','HP','G96','PDih','ImPDih','LJ','Col','Pot','Col-Pro-Pro','LJ-Pro-Pro','Col-Pro-W','LJ-Pro-W','Col-W-W','LJ-W-W']

f_size=15
leg_size=16
tick_size=14

rwt_bool = False

bound_mask = (master_data['d1_mda'] < 12.0) & (master_data['d2_mda'] < 12.0)
cry_mask = (master_data['RMSD-B'] < 2.0)
unbound_mask = ((master_data['d1_mda'] > 12.0) | (master_data['d2_mda'] > 12.0)) & ((master_data['d1_mda'] < 19.0) & (master_data['d2_mda'] < 19.0))

total_bd_data = master_data[bound_mask]
cry_data = master_data[cry_mask]
unb_data = master_data[unbound_mask]


#If want to display Crystal instead of total bound, uncomment this line - 
# bd_data = cry_data
bd_data = total_bd_data

fig,ax = plt.subplots(figsize=(6,4))

ax.plot(bd_data['d1_mda'],bd_data['d2_mda'],label='Bound',linewidth=1.2,linestyle='',color='steelblue',marker='o',markersize=5)
ax.plot(unb_data['d1_mda'],unb_data['d2_mda'],label='Unbound',linewidth=1.2,linestyle='',color='gold',marker='o',markersize=5)

ax.set_xlabel('d1 (nm)',fontsize=f_size)
ax.set_ylabel('d2 (nm)',fontsize=f_size)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
plt.show()

ener_col = ['Bond','HP','G96','PDih','ImPDih','Col-Pro-Pro','LJ-Pro-Pro','Col-Pro-W','LJ-Pro-W','Col-W-W','LJ-W-W']
bound_energy = bd_data[ener_col].sum(axis=1)
unbound_energy = unb_data[ener_col].sum(axis=1)

bound_dict = {}
unbound_dict = {}

if rwt_bool==False:
    bd_data['rbias'][:] = 1.0
    unb_data['rbias'][:] = 1.0

print(bd_data['rbias'])
fig,ax = plt.subplots(figsize=(6,4))
bonded_terms = ['Bond','HP','G96','PDih','ImPDih']
bonded_sum_bd = bd_data[bonded_terms].sum(axis=1)
bonded_sum_unbd = unb_data[bonded_terms].sum(axis=1)
bound_dict['Bonded'] = [np.mean(bonded_sum_bd), np.std(bonded_sum_bd)]
unbound_dict['Bonded'] = [np.mean(bonded_sum_unbd),np.std(bonded_sum_unbd)]
ax.hist(bonded_sum_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
ax.hist(bonded_sum_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3))
ax.set_title('Bonded Energy Distribution',fontsize=f_size)
fig.tight_layout()
plt.show()

fig,ax = plt.subplots(figsize=(6,4))
col_terms = ['Col-Pro-W']
colmb_sum_bd = bd_data[col_terms].sum(axis=1)
colmb_sum_unbd = unb_data[col_terms].sum(axis=1)
bound_dict[col_terms[0]] = [np.mean(colmb_sum_bd), np.std(colmb_sum_bd)]
unbound_dict[col_terms[0]] = [np.mean(colmb_sum_unbd),np.std(colmb_sum_unbd)]

ax.hist(colmb_sum_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
ax.hist(colmb_sum_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.set_title('Col (Pro-W) Energy Distribution',fontsize=f_size)
fig.tight_layout()

fig,ax = plt.subplots(figsize=(6,4))
col_terms = ['Col-Pro-Pro']
colmb_sum_bd = bd_data[col_terms].sum(axis=1)
colmb_sum_unbd = unb_data[col_terms].sum(axis=1)
bound_dict[col_terms[0]] = [np.mean(colmb_sum_bd), np.std(colmb_sum_bd)]
unbound_dict[col_terms[0]] = [np.mean(colmb_sum_unbd),np.std(colmb_sum_unbd)]

ax.hist(colmb_sum_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
ax.hist(colmb_sum_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.set_title('Col (Pro-Pro) Energy Distribution',fontsize=f_size)
fig.tight_layout()
plt.show()

fig,ax = plt.subplots(figsize=(6,4))
lj_terms = ['LJ-Pro-W']
lj_sum_bd = bd_data[lj_terms].sum(axis=1)
lj_sum_unbd = unb_data[lj_terms].sum(axis=1)
bound_dict[lj_terms[0]] = [np.mean(lj_sum_bd), np.std(lj_sum_bd)]
unbound_dict[lj_terms[0]] = [np.mean(lj_sum_unbd),np.std(lj_sum_unbd)]

ax.hist(lj_sum_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
ax.hist(lj_sum_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.set_title('LJ (Pro-W) Energy Distribution',fontsize=f_size)
fig.tight_layout()

fig,ax = plt.subplots(figsize=(6,4))
lj_terms = ['LJ-Pro-Pro']
lj_sum_bd = bd_data[lj_terms].sum(axis=1)
lj_sum_unbd = unb_data[lj_terms].sum(axis=1)
bound_dict[lj_terms[0]] = [np.mean(lj_sum_bd), np.std(lj_sum_bd)]
unbound_dict[lj_terms[0]] = [np.mean(lj_sum_unbd),np.std(lj_sum_unbd)]

ax.hist(lj_sum_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
ax.hist(lj_sum_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.set_title('LJ (Pro-Pro) Energy Distribution',fontsize=f_size)
fig.tight_layout()
plt.show()

fig,ax = plt.subplots(figsize=(6,4))
lj_terms='LJ-W-W'
ljW_bd= bd_data[lj_terms]
ljW_unbd = unb_data[lj_terms]
bound_dict[lj_terms] = [np.mean(ljW_bd), np.std(ljW_bd)]
unbound_dict[lj_terms] = [np.mean(ljW_unbd),np.std(ljW_unbd)]

# ax.hist(ljW_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
# ax.hist(ljW_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.hist(ljW_bd,bins=50,alpha=0.8,label='Bound',color='steelblue',density=True)
ax.hist(ljW_unbd,bins=50,alpha=0.8,label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3),useOffset=False)
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.set_title('LJ Energy (W-W) Distribution',fontsize=f_size)
fig.tight_layout()
plt.show()

fig,ax = plt.subplots(figsize=(6,4))
col_terms = 'Col-W-W'
colW_bd= bd_data[col_terms]
colW_unbd = unb_data[col_terms]
bound_dict[col_terms] = [np.mean(colW_bd), np.std(colW_bd)]
unbound_dict[col_terms] = [np.mean(colW_unbd),np.std(colW_unbd)]

ax.hist(colW_bd,bins=50,alpha=0.8,weights=bd_data['rbias'],label='Bound',color='steelblue',density=True)
ax.hist(colW_unbd,bins=50,alpha=0.8,weights=unb_data['rbias'],label='Unbound',color='gold',density=True)
ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
ax.tick_params(axis='both',labelsize=tick_size)
ax.set_xlabel('Energy (kJ/mol)',fontsize=f_size)
ax.set_ylabel('Frequency',fontsize=f_size)
ax.set_title('Columb Energy (W-W) Distribution',fontsize=f_size)
fig.tight_layout()
plt.show()


output_file = "Bound_EnergyValues.dat"

with open(output_file,'w') as fl:
    fl.write("%-10s\t%-10s\t%-10s\t\t%-10s\t%-10s\n" %('#Energy','Bound-Avg','Bound-STD','Unbound-Avg','Unbound-STD'))
    for term,values in bound_dict.items():
        fl.write("%-10s\t%-10.3f\t%-10.3f\t\t%-10.3f\t%-10.3f\n" %(term,values[0],values[1],unbound_dict[term][0],unbound_dict[term][1]))




# fig,ax = plt.subplots(figsize=(6,4))
# ax.plot(ener_data['Timestep'],ener_data['Bond'],label='Bond',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['HP'],label='HP',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['G96'],label='G96',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['RA'],label='RA',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['PDih'],label='PDih',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['ImPDih'],label='ImPDih',linewidth=1.2)

# ax.set_xlabel('Timestep',fontsize=f_size)
# ax.set_ylabel('Energy (kJ/mol)',fontsize=f_size)
# ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
# ax.tick_params(axis='both',labelsize=tick_size)


# fig,ax = plt.subplots(figsize=(6,4))
# ax.plot(ener_data['Timestep'],ener_data['Col-Pro-Pro'],label='Col: Pro-Pro',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['LJ-Pro-Pro'],label='LJ: Pro-Pro',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['Col-Pro-W'],label='Col: Pro-W',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['LJ-Pro-W'],label='LJ: Pro-W',linewidth=1.2)
# # ax.plot(ener_data['Timestep'],ener_data['LJ-W-W'],label='LJ: W-W',linewidth=1.2)
# ax.set_xlabel('Timestep',fontsize=f_size)
# ax.set_ylabel('Energy (kJ/mol)',fontsize=f_size)
# ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
# ax.tick_params(axis='both',labelsize=tick_size)

# fig,ax=plt.subplots(figsize=(6,4))
# ax.plot(ener_data['Timestep'],ener_data['Col-W-W'],label='Col: W-W',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['LJ-W-W'],label='LJ: W-W',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['Col'],label='Col: System',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['LJ'],label='LJ: System',linewidth=1.2)
# ax.set_xlabel('Timestep',fontsize=f_size)
# ax.set_ylabel('Energy (kJ/mol)',fontsize=f_size)
# ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
# ax.tick_params(axis='both',labelsize=tick_size)


# print("Energy Stats: ")
# print("%-15s%-15s%-15s" %( "Quantity","Mean","Std"))
# for i in range(1,len(col_names)):
#     print("%-15s" %(col_names[i]),end='')
#     print("%-13.2f" %(np.mean(ener_data[col_names[i]])),end='')
#     print("%-13.2f" %(np.std(ener_data[col_names[i]])),end='')
#     print("\n")


# subset = ['Bond','HP','G96','RA','PDih','ImPDih','Col-Pro-Pro','LJ-Pro-Pro','Col-Pro-W','LJ-Pro-W','Col-W-W','LJ-W-W']
# total_sum = ener_data[subset].sum(axis=1)

# fig,ax = plt.subplots(figsize=(6,4))
# ax.plot(ener_data['Timestep'],total_sum,label='Summed Energy',linewidth=1.2)
# ax.plot(ener_data['Timestep'],ener_data['Pot'],label='Potential Energy',linewidth=1.2)
# ax.set_xlabel('Timestep',fontsize=f_size)
# ax.set_ylabel('Energy (kJ/mol)',fontsize=f_size)
# ax.legend(fontsize=leg_size,fancybox=True,framealpha=0.6,markerscale=1.2)
# ax.tick_params(axis='both',labelsize=tick_size)
# plt.show()