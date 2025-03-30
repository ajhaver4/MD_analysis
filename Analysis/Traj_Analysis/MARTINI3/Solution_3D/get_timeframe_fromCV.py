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

CV_flag = sys.argv[2]

CV_val = sys.argv[3]
output_lbl = sys.argv[4]

tolerance = 0.001

mask = (dist_data[CV_flag] >= float(CV_val) - tolerance) & (dist_data[CV_flag] <= float(CV_val) + tolerance)

timeframe = dist_data[mask]['Timestep'].values
cv_values = dist_data[mask][CV_flag].values

print(len(timeframe))
print(cv_values)


#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
outfile = output_lbl + ".pickle"
pick_path = outfile

#Storing the analysis in pickle file
with open(pick_path,"wb") as pick_handle:
    pickle.dump(timeframe,pick_handle)