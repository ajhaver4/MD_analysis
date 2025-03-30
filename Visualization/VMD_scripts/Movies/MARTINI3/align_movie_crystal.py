import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances as dist
from MDAnalysis.lib.distances import distance_array as d_array
import math
import MDAnalysis.analysis.rms as rms
import time as time_mod
from MDAnalysis.analysis import align
from MDAnalysis.core.universe import Merge


template_univ = mda.Universe(sys.argv[1],sys.argv[2])
mobile_univ = mda.Universe(sys.argv[3],sys.argv[4])

print(template_univ.atoms.positions.shape)
print(mobile_univ.atoms.positions.shape)
aligner_B = align.AlignTraj(mobile_univ, template_univ,select ='bynum 1-469',in_memory=True).run()    

file_str = "Aligned" + sys.argv[4]

#Check if trajectory 
if sys.argv[4].endswith('.xtc'):
    mobile_univ.atoms.write(file_str, frames='all')
else:
    mobile_univ.atoms.write(file_str)