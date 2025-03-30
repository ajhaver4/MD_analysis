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

print("No. of frames: ",n_frames)
print("Current time: ", u.trajectory.time)
print("Length of trajectory (us): ",u.trajectory.totaltime)

"""
    --------------------------------------- RMSF -----------------------------------
    """

avg_struct = align.AverageStructure(u,u,select='bynum 1-916',ref_frame=0).run()
ref_struct = avg_struct.universe 
ref_group = ref_struct.atoms
ref_group.write(outfile_pref+"_ChainA_average.pdb")
aligner_A = align.AlignTraj(u, ref_struct,select='bynum 1-458 and name BB',in_memory=True).run()

rmsf_chainA_sel = u.select_atoms('bynum 1-458 and name BB')
# rmsf_chainB_sel = u.select_atoms('bynum 470-938 and name BB')
rmsf_chainA = rms.RMSF(rmsf_chainA_sel).run()
# rmsf_chainBtoA = rms.RMSF(rmsf_chainB_sel).run()



avg_structB = align.AverageStructure(u_copy,u_copy,select='bynum 1-916',ref_frame=0).run()
ref_structB = avg_structB.universe
ref_broupB = ref_structB.atoms
ref_broupB.write(outfile_pref+"_ChainB_average.pdb")
print(ref_structB.atoms)
aligner_B = align.AlignTraj(u_copy, ref_structB,select='bynum 459-916 and name BB',in_memory=True).run()    #The ref structure for B only contains 469 atoms of chain B. 

# rmsf_chainA_sel = u_copy.select_atoms('bynum 1-469 and name BB')
rmsf_chainB_sel = u_copy.select_atoms('bynum 459-916 and name BB')
rmsf_chainB = rms.RMSF(rmsf_chainB_sel).run()
# rmsf_chainAtoB = rms.RMSF(rmsf_chainA_sel).run()


outfileA = outfile_pref + "_RMSF_chainA.dat"
with open(outfileA, 'a') as fl3:
    fl3.write("#Residue\tRMSF\n")
    for i in range(len(rmsf_chainA_sel.atoms)):
        fl3.write("%d\t%.3f" %(i,rmsf_chainA.rmsf[i]))
        fl3.write("\n")

        

    fl3.write("\n")

outfileB = outfile_pref + "_RMSF_chainB.dat"
with open(outfileB, 'a') as fl3:
    fl3.write("#Residue\tRMSF\n")
    
    for i in range(len(rmsf_chainB_sel.atoms)):
        fl3.write("%d\t%.3f" %(i,rmsf_chainB.rmsf[i]))
        fl3.write("\n")

    fl3.write("\n")
