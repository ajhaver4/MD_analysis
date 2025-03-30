import sys
import numpy as np
import pandas
from scipy import constants


"""
Usage
python sep_hills_blocksMem.py HILLS_file history_bool block_start time_int

HILLS_file: HILLS file to separate in blocks
history_bool: True or False. If True, the script will keep the history of the previous blocks
block_start: Percentage of the total time to start splitting the blocks
time_int: Time interval for the blocks in microseconds
"""

data = pandas.read_csv(sys.argv[1],comment='#',names=['time','d1','d2','sigma1','sigma2','height','biasf'],delimiter='\s+')

print(data)
history_bool = sys.argv[2]
#Dividing into blocks
block_start = float(sys.argv[3])     #Percentage of the total time

init_time = 0
# time_points = [90,236]
time_int = float(sys.argv[4]) #us


header_col = []
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('#'):
            header_col.append(line)
        else:
            break


raw_time = len(data['time'])*300*1e-6
print("Raw time of HILLS file:", raw_time)


timesteps_array = np.arange(0,len(data['time']),1)
tau=300
#Time points to exclude from the analysis
mask1 = (timesteps_array > (33.0*1e6/tau)) & (timesteps_array < 66.0*1e6/tau)
mask2 = (timesteps_array > 95.2*1e6/tau) & (timesteps_array < 97.8*1e6/tau)
mask3 = (timesteps_array > 101.95*1e6/tau) & (timesteps_array < 103.3*1e6/tau)
mask4 = (timesteps_array > 127.6*1e6/tau) & (timesteps_array < 128.5*1e6/tau)
mask5 = (timesteps_array > 139*1e6/tau) & (timesteps_array < 142.5*1e6/tau)
mask = mask1 | mask2 | mask3 | mask4 | mask5

real_data = data[~mask]
total_time = len(real_data['time'])*300*1e-6
time_points = np.arange(total_time*block_start,total_time,time_int)
print("Total time of HILLS file:", total_time)
#Check the difference between the last time point adn the total time
diff = total_time - time_points[-1]
#If the difference is less than 50% of the block size, add the remaining steps in the last block 
if diff < 0.5*time_int:
    time_points[-1] = total_time
else:
    time_points = np.append(time_points,total_time)

print("Time points for blocks:", time_points)

timesteps_array = np.arange(0,len(real_data['time']),1)

tau=300
step_min = int(init_time*1e6/tau)
for i in range(len(time_points)):
    
    step_max = int(time_points[i]*1e6/tau)
    time_label = int(time_points[i])
    # if i == len(time_points)-1:
    #     step_max = len(real_data)
    #     time_label = int(step_max*300*1e-6)
    
    filter_mask = (timesteps_array >= step_min) & (timesteps_array <= step_max)
    # final_mask = ~mask & filter_mask
    final_mask = filter_mask 
    print("Block:", i)
    # block_data = data[final_mask]
    block_data = real_data[final_mask]
    
    # block_data.to_csv(f"HILL_block_{time_label}us.dat",index=False,sep='\t')

    with open(f"HILL_block_{time_label}us.dat", 'w') as f:
        for line in header_col:
            f.write(line)
        block_data.to_csv(f,sep='\t',index=False,header=False)

    if history_bool==False:
        print("Discarding previous history")
        step_min = step_max


