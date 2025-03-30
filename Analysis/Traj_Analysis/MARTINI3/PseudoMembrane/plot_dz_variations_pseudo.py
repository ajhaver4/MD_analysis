import numpy as np 
import pandas
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.collections import LineCollection

file=sys.argv[1]
dz_data = pandas.read_csv(file,comment='#',sep='\t')
master_data = pandas.read_csv(sys.argv[2],comment='#',sep='\t')


#Remove duplicates 
dist_data = dz_data.drop_duplicates(subset=['Timestep'],ignore_index=True)
ener_data = master_data.drop_duplicates(subset=['Timestep'],ignore_index=True)

data = dz_data.join(master_data.set_index('Timestep'),on='Timestep',how='inner',lsuffix = '_mda',rsuffix = '_plumed')


# disp_ts = 100   #ns
# interval = int(disp_ts*1e3 / 200)  #200 ps time interval between data points in the data file
filter_data = data 


bd_mask = (filter_data['d1'] < 13.0) & (filter_data['d2'] < 13.0)
# cry_mask = filter_data['RMSD-BtoA'] < 2.0
unbound_mask = ((filter_data['d1'] > 13.0) | (filter_data['d2'] > 13.0)) & ((filter_data['d1'] < 19.0) & (filter_data['d2'] < 19.0))

bd_data = filter_data[bd_mask]
unb_data = filter_data[unbound_mask]

fig,[ax,ax1] = plt.subplots(2,1)

ax.hist(bd_data['COMz_chainA'],bins=50,label=r'$z_{disp}|Bd$')
# ax.hist(bd_data['z_A_membres'],bins=50,label=r'$z_{surf}|Bd$')
ax.set_xlabel("z- distance")
ax.set_ylabel("P(z)")
ax.legend()

ax1.hist(unb_data['COMz_chainB'],bins=50,label=r'$z_{disp}|Unb$')
# ax1.hist(unb_data['z_A_membres'],bins=50,label=r'$z_{surf}|Unb$')
ax1.set_xlabel("z- distance")
ax1.set_ylabel("P(z)")
ax1.legend()
fig.suptitle(r"$z_{disp} ~vs~ z_{surf}$")

plt.show()
def integrate_deltaz(xvalues,filter=False):
    if filter:
        mask = (xvalues >0.5) & (xvalues < 2.0)
        xhisto = np.histogram(xvalues[mask],bins=100)
    else:
        xhisto = np.histogram(xvalues,bins=100)
    dz = xhisto[1][1]-xhisto[1][0]    #bin size
    integral = np.sum(dz*xhisto[0])/np.max(xhisto[0])    # Sum of all values in the integration area = Ni*bin_size
    return(integral)



#Integral over chain A


bd_dz_dispA = integrate_deltaz(bd_data['COMz_chainA'])
bd_dz_dispB = integrate_deltaz(bd_data['COMz_chainB'])

unb_dz_dispA = integrate_deltaz(unb_data['COMz_chainA'])
unb_dz_dispB = integrate_deltaz(unb_data['COMz_chainB'])


print("dz integrals")

print("Bound: ")
print(bd_dz_dispA)
print(bd_dz_dispB)
print("Unbound: ")
print(unb_dz_dispA)
print(unb_dz_dispB)