import sys
import numpy as np
import pandas
from scipy import constants


data = pandas.read_csv(sys.argv[1],comment='#',names=['time','d1','d2','sigma1','sigma2','height','biasf'],delimiter='\s+')

history_bool = sys.argv[2]
tau = int(sys.argv[3])
time_int = int(sys.argv[4])
wt_min = int(sys.argv[5])
wt_max = int(sys.argv[6])

#Dividing into blocks

def parse_input(file):
    with open(file, 'r') as fl:
        lines = f.readlines()

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

wt_flag = False
step_min = int(init_time*1e6/tau)
for i in range(len(time_points)):
    
    step_max = int(time_points[i]*1e6/tau)
    time_label = int(time_points[i])

    print("Block:", i)

    

    if time_points[i]> wt_min:
        step_max = int(wt_min*1e6/tau)
        wt_step_min = step_max
        wt_flag = True

        block_data = data[step_min:step_max]

        # block_data.to_csv(f"HILL_block_{time_points[i]}us.dat",index=False,sep='\t')
        print("Writing block data to HILLS:", time_label)
        with open(f"HILL_block_{time_label}us.dat", 'w') as f:
            for line in header_col:
                f.write(line)
            block_data.to_csv(f,sep='\t',index=False,header=False)

        
        if time_points[i]< wt_max:
            wt_step_max = int(time_points[i]*1e6/tau)

            wt_block_data = data[wt_step_min:wt_step_max]
            print("Writing WT Hills to file from time %f - %f" %(wt_step_min*tau/1e6, wt_step_max*tau/1e6))
            with open(f"HILL_block_{time_label}us_WT.dat", 'w') as f:
                for line in header_col:
                    f.write(line)
                wt_block_data.to_csv(f,sep='\t',index=False,header=False)
            

        elif time_points[i] >= wt_max:
            wt_flag = False
            wt_step_max = int(wt_max*1e6/tau)
            final_step_max = int(time_points[i]*1e6/tau)

            wt_block_data = data[wt_step_min:wt_step_max]

            last_block_data = data[wt_step_max:final_step_max]
            std_block_data = pandas.concat([block_data,last_block_data])

            # print(block_data)
            # print(last_block_data)
            # print(std_block_data)
            print("Writing WT Hills to file from time %f - %f" %(wt_step_min*tau/1e6, wt_step_max*tau/1e6))
            with open(f"HILL_block_{time_label}us_WT.dat", 'w') as f:
                for line in header_col:
                    f.write(line)
                wt_block_data.to_csv(f,sep='\t',index=False,header=False)
            print("Writing standard Hills to file from time %f - %f" %(wt_step_max*tau/1e6, final_step_max*tau/1e6))
            with open(f"HILL_block_{time_label}us.dat", 'w') as f:
                for line in header_col:
                    f.write(line)
                std_block_data.to_csv(f,sep='\t',index=False,header=False)
            

    else:       
  
    
        block_data = data[step_min:step_max]
        # block_data.to_csv(f"HILL_block_{time_points[i]}us.dat",index=False,sep='\t')
        with open(f"HILL_block_{time_label}us.dat", 'w') as f:
            for line in header_col:
                f.write(line)
            block_data.to_csv(f,sep='\t',index=False,header=False)


    if history_bool==False:
        print("Discarding previous history")
        step_min = step_max


