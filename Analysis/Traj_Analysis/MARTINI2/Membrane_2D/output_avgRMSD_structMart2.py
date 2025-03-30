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


#Reading in all the inputs
struct_file = sys.argv[1]
traj_file = sys.argv[2]
outfile_pref = sys.argv[3]

#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()

#Selecting only protein beads
protein = u.select_atoms("byres name BB")
protein_rms = u.select_atoms("byres name BB")
protein_rmsf = u.select_atoms("byres name BB")
print("Total mass: ",protein.total_mass())
# sys.exit()
#print(protein)

#Selecting protein chains
pro_A = protein.select_atoms("bynum 1-458")
pro_B = protein.select_atoms("bynum 459-916")
pro_B_mic = protein_rms.select_atoms("bynum 459-916")
atom_len = len(pro_A)
pro_A_bb = pro_A.select_atoms("name BB")
pro_B_bb = pro_B.select_atoms("name BB")

n_frames = u.trajectory.n_frames
ts = u.trajectory.ts
print("Timestep:",ts.dt)
global box_dims
box_dims = ts.dimensions
print("Box dimension: ", box_dims)
print("Number of frames: ", n_frames)



#First get average coordinates of chain A by aligning all frames to the first frame
#RMS-aligns all trajectory frames to the first frame -> The alignment removes translation and rotation motions of the protein to looks at the internal conformations 
#Chain A is selected bynum 1-469
avg_struct = align.AverageStructure(u,u,select='bynum 1-458',ref_frame=0).run()
ref_struct = avg_struct.universe 
ref_group = ref_struct.atoms

#Now do the same for chain B
#Use a different universe since the previous universe has been modified
avg_structB = align.AverageStructure(u_copy,u_copy,select='bynum 459-916',ref_frame=0).run()
ref_structB = avg_structB.universe
ref_broupB = ref_structB.atoms

#Merge the averagted coordinates the two chains and output to a pdb file which will be then used to calculate the RMSD
merge_ref = Merge(ref_struct.universe.atoms,ref_structB.universe.atoms)
merge_ref.atoms.write(outfile_pref+"Averaged_Struct_chainAB.pdb")