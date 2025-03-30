import sys
import numpy as np
import pandas
from scipy import constants


data = pandas.read_csv(sys.argv[1],comment='#',names=['time','d1','d2','sigma1','sigma2','height','biasf'],delimiter='\s+')

history_bool = sys.argv[2]
tau = int(sys.argv[3])
time_int = int(sys.argv[4])
#Dividing into blocks

header_col = []
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('#'):
            header_col.append(line)
        else:
            break

init_time = 0
# time_points = [30,73,106,163,208]


total_time = len(data['time'])*tau*1e-6
# time_int = 20 #us
time_points = np.arange(total_time*0.5,total_time,time_int)
print("Total time of HILLS file:", total_time)

#Check the difference between the last time point adn the total time
diff = total_time - time_points[-1]
#If the difference is less than 50% of the block size, add the remaining steps in the last block 
if diff < 0.5*time_int:
    time_points[-1] = total_time
else:
    time_points = np.append(time_points,total_time)

print("Time points for blocks:", time_points)


step_min = int(init_time*1e6/tau)
for i in range(len(time_points)):
    
    step_max = int(time_points[i]*1e6/tau)
    time_label = int(time_points[i])
    # if i == len(time_points)-1:
    #     #Check how many steps are left after the last block
    #     last_steps = len(data)

    #     diff = (last_steps - step_max) * 300 * 1e-6

    #     #If the difference is less than 50% of the block size, add the remaining steps to the last block
    #     if diff < 0.5*time_int:
    #         step_max = len(data)
    #         time_label = int(step_max*300*1e-6)
    #     else:
    #         #If the difference is greater than 50% of the block size, store the new block 
    #         store_last_bool = True
    #     step_max = len(data)
    #     time_label = int(step_max*300*1e-6)
    

    print("Block:", i)
    block_data = data[step_min:step_max]
    # block_data.to_csv(f"HILL_block_{time_points[i]}us.dat",index=False,sep='\t')
    with open(f"HILL_block_{time_label}us.dat", 'w') as f:
        for line in header_col:
            f.write(line)
        block_data.to_csv(f,sep='\t',index=False,header=False)

    if history_bool==False:
        print("Discarding previous history")
        step_min = step_max


