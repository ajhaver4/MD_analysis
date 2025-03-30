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
wall_max = float(sys.argv[4])


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
for i in range(len(states_data['st'])):

    coords = [float(states_data['xmin'][i]),float(states_data['xmax'][i]),float(states_data['ymin'][i]),float(states_data['ymax'][i])]
    states[states_data['st'][i]] = coords
    dG[states_data['st'][i]] = []

#Create a grid of coordinates
d1_grid,d2_grid = np.meshgrid(d1_bins,d2_bins)
#Create a 2D array for biaspotential
bias_pot = np.zeros((len(d1_bins),len(d2_bins)))





def integrate_fes(bias_pot,states,dG,d1_bins,d2_bins):
    k = constants.value('Boltzmann constant')
    Ava_no = constants.value('Avogadro constant')
    temp = 310
    kbt_ = 1000/(k*temp*Ava_no)   #kJ/mol

    dx = d1_bins[1]-d1_bins[0]
    dy = d2_bins[1]-d2_bins[0]
    Z_states = {}
    for st,coords in states.items():
        xmin=coords[0]
        if xmin==-1.0:
            xmin=np.min(d1_bins)+0.5
        xmax=coords[1]
        if xmax==-1.0:
            xmax=wall_max
        ymin=coords[2]
        if ymin==-1.0:
            ymin=np.min(d2_bins)+0.5
        ymax=coords[3]
        if ymax==-1.0:
            ymax=wall_max

        mask1 = (d1_bins>=xmin) & (d1_bins<=xmax)
        mask2 = (d2_bins>=ymin) & (d2_bins<=ymax)

        bias1= bias_pot[mask1,:]
        bias2=bias1[:,mask2]

        Z = np.sum(np.exp(-1*bias2/kbt_)*dx*dy)
        Z_states[st]=Z

    #Calculating dG
    dG_state = 0
    Z_unb = Z_states['unb']
    for st,Z_bd in Z_states.items():
        dG_state = -kbt_*np.log(Z_bd/Z_unb)
        dG[st].append(dG_state)

    return(dG)



print("Reading HILLS file and building FES")

timesteps = []
#Averaging time
tau=280   #Deposition time
t_max = len(data['time'])*tau
print("Max time: ",t_max)
t_avg = float(sys.argv[3])*1e6

print("Averagin time start..",t_avg)

for idx in range(len(data['time'])):
    time = idx*tau
    curr_d1=data['d1'][idx]
    curr_d2=data['d2'][idx]
    h = data['height'][idx]
    s1 = data['sigma1'][idx]
    s2 = data['sigma2'][idx]

    if time > t_avg:
        avg_scaling = (t_max-time)/(t_max-t_avg)
    else:
        avg_scaling=1


    # for i in range(len(d1_grid)-1):
    #     # print("Bin no.: ",i)
    #     for j in range(len(d2_grid)-1):
    #         x=(d1_grid[i,j]+d1_grid[i+1,j])/2   #Bin center
    #         y=(d2_grid[i,j]+d2_grid[i,j+1])/2   #Bin center
    #         bias_x = (x-curr_d1)**2/(2*(s1**2))
    #         bias_y = (y-curr_d2)**2/(2*(s2**2))
    #
    #         bias_pot[i,j]+= -1*h*np.exp(-1*(bias_x+bias_y))
            # print(bias_pot[i,j])

    bias_x = ((d1_grid+bin_width_d1*0.5) - curr_d1)**2/(2*(s1**2))
    bias_y = ((d2_grid+bin_width_d2*0.5) - curr_d2)**2/(2*(s2**2))

    expo = np.exp(-1*(bias_x+bias_y))


    bias_pot+= -1*h*avg_scaling*expo        #Minus sign

    dG = integrate_fes(bias_pot,states,dG,d1_bins,d2_bins)
    timesteps.append(time)

    # if data['time'][idx]>100000:
    #     break

    if idx%100000==0:
        print("Current time: %.2f us" %(time/1e6))
    # sys.exit()

print("Finished building FES...")
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
