import sys
import numpy as np
import pandas
from scipy import constants

"""
Sample code to just check meshgrid format and inddexing
"""
# x = np.linspace(0,10,11)
# y= np.linspace(-10,0,11)
#
# x_grid,y_grid = np.meshgrid(x,y,indexing='ij')
#
#
# print(x_grid)
# print(y_grid)
#
# for i in range(len(x_grid)):
#     for j in range(len(y_grid)):
#         print(x_grid[i,j],y_grid[i,j])

"""
Mian code
"""
data = pandas.read_csv(sys.argv[1],comment='#',names=['time','d1','d2','sigma1','sigma2','height','biasf'],delimiter='\s+')
states_data = pandas.read_csv(sys.argv[2],delimiter='\s+')

avg_bool=True
if float(sys.argv[3]) == -1.0:
    avg_bool=False
    t_avg=-1
else:
    t_avg = float(sys.argv[3])*1e6



k = constants.value('Boltzmann constant')
Ava_no = constants.value('Avogadro constant')
temp = 310
kbt = (k*temp*Ava_no)/1000   #kJ/mol

print("Vakue of kT : %.3f kJ/mol" %kbt)

#Create a grid

#Get number of bins for d1 and d2
min_d1 = np.amin(data['d1'])-0.5
min_d2 = np.amin(data['d2'])-0.5
max_d1 = np.amax(data['d1'])+1.0
max_d2 = np.amax(data['d2'])+1.0

d1_bins = np.linspace(min_d1,max_d1,800)
d2_bins = np.linspace(min_d2,max_d2,800)

bin_width_d1 = d1_bins[1]-d1_bins[0]
bin_width_d2 = d2_bins[1]-d2_bins[0]

#Creating states dictionary to hold boudnaries
states = {}
dG = {}

dG_unb_calc={}
unbound_labels = []
for i in range(len(states_data['st'])):

    coords = [float(states_data['xmin'][i]),float(states_data['xmax'][i]),float(states_data['ymin'][i]),float(states_data['ymax'][i])]
    states[states_data['st'][i]] = coords
    dG[states_data['st'][i]] = []


    if 'unb' in states_data['st'][i]:
        unbound_labels.append(states_data['st'][i])

#Create a grid of coordinates
d1_grid,d2_grid = np.meshgrid(d1_bins,d2_bins)
#Create a 2D array for biaspotential
bias_pot = np.zeros((len(d1_bins),len(d2_bins)))




def integrate_fes(fesgrid,states,dG,dG_unb_calc,d1_bins,d2_bins):

    dx = d1_bins[1]-d1_bins[0]
    dy = d2_bins[1]-d2_bins[0]
    Z_states = {}
    for st,coords in states.items():
        xmin=coords[0]
        if xmin==-1.0:
            xmin=np.min(d1_bins)+0.5
        xmax=coords[1]
        if xmax==-1.0:
            xmax=np.max(d1_bins)-1.0
        ymin=coords[2]
        if ymin==-1.0:
            ymin=np.min(d2_bins)+0.5
        ymax=coords[3]
        if ymax==-1.0:
            ymax=np.max(d2_bins)-1.0

        mask1 = (d1_bins>=xmin) & (d1_bins<=xmax)
        mask2 = (d2_bins>=ymin) & (d2_bins<=ymax)

        if 'unb' in st:
            F1 = fesgrid[mask2,:]
            F2 = fesgrid[:,mask1]

            F12 = F1[:,mask1]

            Z = np.sum(np.exp(-1*F1/kbt)*dx*dy) + np.sum(np.exp(-1*F2/kbt)*dx*dy) - np.sum(np.exp(-1*F12/kbt)*dx*dy)
        
        else:
            F= fesgrid[mask2,:]
            F_of_s=F[:,mask1]
            Z = np.sum(np.exp(-1*F_of_s/kbt)*dx*dy)

        # F= fesgrid[mask1,:]
        # F_of_s=F[:,mask2]

        # Z = np.sum(np.exp(-1*F_of_s/kbt)*dx*dy)
        Z_states[st]=Z

    #Calculating dG
    dG_state = 0
    Z_unb = Z_states['unb']

    
    for st,Z_bd in Z_states.items():
        dG_state = -kbt*np.log(Z_bd/Z_unb)
        dG[st].append(dG_state)

        # for unb_st in unbound_labels:
        #     dG_state_unb = -kbt*np.log(Z_bd/Z_states[unb_st])
        #     lbl = st+"_"+unb_st
        #     if lbl not in dG_unb_calc.keys():
        #         dG_unb_calc[lbl] = []
            
        #     dG_unb_calc[lbl].append(dG_state_unb)
        

    return(dG)


def integrate_1D(fes_grid,d1_bins,d2_bins):

    d1_fes = np.zeros(len(d1_bins))
    d2_fes = np.zeros(len(d2_bins))

    d1_bias = np.zeros(len(d1_bins))
    d2_bias = np.zeros(len(d2_bins))

    dx = d1_bins[1]-d1_bins[0]
    dy = d2_bins[1]-d2_bins[0]

    for i in range(len(d1_bins)):
        d1_fes[i] = -kbt*np.log(np.sum(np.exp(-1*fes_grid[i,:]/kbt)*dy))
        d1_bias[i] = np.sum(fes_grid[i,:])

    for i in range(len(d2_bins)):
        d2_fes[i] = -kbt*np.log(np.sum(np.exp(-1*fes_grid[:,i]/kbt)*dx))
        d2_bias[i] = np.sum(fes_grid[:,i])

    return(d1_fes,d2_fes,d1_bias,d2_bias)

print("Reading HILLS file and building FES")

timesteps = []
#Averaging time
tau=300   #Deposition time
t_max = len(data['time'])*tau
print("Max time: ",t_max)


print("Averagin time start..",t_avg)
avg_scaling=1

for idx in range(len(data['time'])):
    time = idx*tau
    curr_d1=data['d1'][idx]
    curr_d2=data['d2'][idx]
    h = data['height'][idx]
    s1 = data['sigma1'][idx]
    s2 = data['sigma2'][idx]

    if time > 33000000 and time < 65500000:
        continue

    if time > 95500000 and time < 98000000:
        continue

    if time > 139000000 and time < 142500000:
        continue

    if avg_bool:
        if time > t_avg:
            avg_scaling = (t_max-time)/(t_max-t_avg)
        else:
            avg_scaling=1

    dx = d1_bins[1]-d1_bins[0]
    dy = d2_bins[1]-d2_bins[0]
    bias_x = ((d1_grid+bin_width_d1*0.5) - curr_d1)**2/(2*(s1**2))
    bias_y = ((d2_grid+bin_width_d2*0.5) - curr_d2)**2/(2*(s2**2))
    expo = np.exp(-1*(bias_x+bias_y))
    
    bias_pot+= -1*h*avg_scaling*expo

    
    dG = integrate_fes(bias_pot,states,dG,dG_unb_calc,d1_bins,d2_bins)
    timesteps.append(time)

    # if data['time'][idx]>100000:
    #     break

    if idx%100000==0:
        print("Current time: %.2f us" %(time/1e6))
      
    # sys.exit()

print("Finished building FES...")
print("Calculating 1D FES...")
d1_fes,d2_fes,d1_bias,d2_bias = integrate_1D(bias_pot,d1_bins,d2_bins)
print("Writing to file...")

filename = "fes_build_" + str(sys.argv[3]) + "_us.dat"
with open(filename,"w") as fl:
    fl.write("#! FIELDS d1 d2 free")
    fl.write("\n")
    for i in range(len(d1_grid)):
        for j in range(len(d2_grid)):
            fl.write("%2.6f\t%2.6f\t%.6f\n" %(d1_grid[i,j],d2_grid[i,j],bias_pot[i,j]))

with open("dG_vs_time.dat","w") as fl1:
    fl1.write("#Timestep")
    for state in dG:
        fl1.write("\t%s" %(state))
    fl1.write("\n")

    for i in range(len(timesteps)):
        fl1.write(str(timesteps[i]))
        for st,deltaG in dG.items():
            fl1.write("\t%11.3f" %(deltaG[i]))
        fl1.write("\n")

# with open("dG_unb_vs_time.dat","w") as fl1:
#     fl1.write("#Timestep")
#     for state in dG_unb_calc:
#         fl1.write("\t%s" %(state))
#     fl1.write("\n")

#     for i in range(len(timesteps)):
#         fl1.write(str(timesteps[i]))
#         for st,deltaG in dG_unb_calc.items():
#             fl1.write("\t%11.3f" %(deltaG[i]))
#         fl1.write("\n")
with open("fes_d1.dat","w") as fl1:
    fl1.write("#d1\tFES\tbias\n")
    for i in range(len(d1_fes)):
        fl1.write(str(d1_bins[i]))
        fl1.write("\t")
        fl1.write(str(d1_fes[i]))
        fl1.write("\t")
        fl1.write(str(d1_bias[i]))
        fl1.write("\n")

with open("fes_d2.dat","w") as fl1:
    fl1.write("#d2\tFES\tbias\n")
    for i in range(len(d2_fes)):
        fl1.write(str(d2_bins[i]))
        fl1.write("\t")
        fl1.write(str(d2_fes[i]))
        fl1.write("\t")
        fl1.write(str(d2_bias[i]))
        fl1.write("\n")


