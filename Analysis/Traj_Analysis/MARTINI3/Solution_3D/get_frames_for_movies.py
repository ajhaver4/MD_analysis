import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.markers import MarkerStyle as markerstyle
import math
from scipy.interpolate import griddata
import pandas
from scipy import constants
import matplotlib.patches as patches
import matplotlib.cm as cm

dist_data = pandas.read_csv(sys.argv[1],comment='#',sep='\t')    #Master distances file

CV_flag = sys.argv[4]
time_range = (int(sys.argv[2])*1e6,int(sys.argv[3])*1e6)

mask1 = (dist_data['Timestep'] >= time_range[0]) & (dist_data['Timestep'] <= time_range[1])

mask_bnd = (dist_data['d1'] <= 13.0 ) & (dist_data['d2']<=13.0)
unb_data = dist_data[~mask_bnd]


dist_data = dist_data[mask1]

print("No. of frames in time range: ",len(dist_data['Timestep']))
dist_data.reset_index(drop=True,inplace=True)
unb_data.reset_index(drop=True,inplace=True)

crystal_times = []
indices = []

for i in range(len(dist_data['Timestep'])):

    if dist_data[CV_flag][i] <= 2.0:
        crystal_times.append(dist_data['Timestep'][i])
        indices.append(i)


n_crystal = len(crystal_times)
print("Total number of frames in crystal: ",n_crystal)
print("Crystal Bound state range: ", crystal_times[0], crystal_times[-1])

movie_times = [dist_data['Timestep'][indices[0]-n_crystal], dist_data['Timestep'][indices[-1]+n_crystal]]
print("Time range for movie : ",movie_times)


unb_range1 = unb_data['Timestep'][unb_data['Timestep'] <= crystal_times[0]].to_numpy()
unb_range2 = unb_data['Timestep'][unb_data['Timestep'] >= crystal_times[-1]].to_numpy()
#Closest unbound frame to crystal structure

unb_min = np.argmin(np.abs(unb_range1 - crystal_times[0]))
unb_max = np.argmin(np.abs(unb_range2 - crystal_times[-1]))

print("Unbound min: ",unb_range1[unb_min])
print("Unbound max: ",unb_range2[unb_max])

print("Time range for movie :", unb_range1[unb_min], unb_range2[unb_max])
print("No. of frames :", (unb_range2[unb_max]-unb_range1[unb_min])/200)