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

rmsd_flag = sys.argv[2]
output_lbl = sys.argv[3]

rmsd1 = 0.6 
rmsd2 = 2.0 


mask1 = (dist_data[rmsd_flag] <= float(rmsd1))
mask2 = (dist_data[rmsd_flag] > float(rmsd1)) & (dist_data[rmsd_flag] <= float(rmsd2))



timeframes1 = dist_data[mask1]['Timestep'].values
timeframes2 = dist_data[mask2]['Timestep'].values

rmsd_values1 = dist_data[mask1][rmsd_flag].values
rmsd_values2 = dist_data[mask2][rmsd_flag].values


print("Frames for RMSD <= 0.6 nm: ", len(timeframes1))
print("Frames for 0.6 < RMSD <= 2.0 nm: ", len(timeframes2))



#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
outfile1 = output_lbl + "_0p6nm.pickle"
outfile2 = output_lbl + "_2nm.pickle"
pick_path1 = outfile1
pick_path2 = outfile2

#Storing the analysis in pickle file
with open(pick_path1,"wb") as pick_handle:
    pickle.dump(timeframes1,pick_handle)

with open(pick_path2,"wb") as pick_handle:
    pickle.dump(timeframes2,pick_handle)