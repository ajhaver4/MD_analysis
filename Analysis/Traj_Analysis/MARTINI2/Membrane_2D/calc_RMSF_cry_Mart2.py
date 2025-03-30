import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
# from MDAnalysis.analysis import distances as dist
from MDAnalysis.lib.distances import distance_array as d_array
import math
import MDAnalysis.analysis.rms as rms
import time as time_mod
from MDAnalysis.analysis import align
from MDAnalysis.core.universe import Merge

import MDAnalysis.analysis.distances as distances

#Reading in all the inputs
struct_file = sys.argv[1]
traj_file = sys.argv[2]

cry_univ = mda.Universe(sys.argv[3],sys.argv[4])
outfile_pref = sys.argv[5]
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

pro_A_cry = cry_univ.select_atoms("bynum 1-458")
pro_B_cry = cry_univ.select_atoms("bynum 459-916")

pro_B_mic = protein_rms.select_atoms("bynum 459-916")
atom_len = len(pro_A)

pro_A_bb = pro_A.select_atoms("name BB")
pro_B_bb = pro_B.select_atoms("name BB")

proA_bb_cry = pro_A_cry.select_atoms("name BB")
proB_bb_cry = pro_B_cry.select_atoms("name BB")

n_frames = u.trajectory.n_frames
ts = u.trajectory.ts
print("Timestep:",ts.dt)
global box_dims
box_dims = ts.dimensions
print("Box dimension: ", box_dims)
print("Number of frames: ", n_frames)

def check_image(pos1,pos2):
    del_x = pos1[0] - pos2[0]
    del_y = pos1[1] - pos2[1]
    del_z = pos1[2] - pos2[2]

    #print("Del_y:", del_y)
    #print("Del_z:", del_z)

    transVec = np.zeros((3))
    mic = False
    if abs(del_x) > box_dims[0]/2:
        mic = True
        #print("X")
        if del_x > 0:
            transVec[0] = box_dims[0]
        else:
            transVec[0] = -box_dims[0]
    if abs(del_y) > box_dims[1]/2:
        mic = True
        #print("Y")
        if del_y > 0:
            transVec[1] = box_dims[1]
        else:
            transVec[1] = -box_dims[1]
    if abs(del_z) > box_dims[2]/2:
        mic = True
        #print("Z")
        if del_z > 0:
            transVec[2] = box_dims[2]
        else:
            transVec[2] = -box_dims[2]

    r_nopbc = ((del_x)**2 + (del_y)**2 + (del_z)**2)**0.5

    del_x = del_x - (box_dims[0])*round(del_x/box_dims[0],0)
    del_y = del_y - (box_dims[1])*round(del_y/box_dims[1],0)
    del_z = del_z - (box_dims[2])*round(del_z/box_dims[2],0)
    r = ((del_x)**2 + (del_y)**2 + (del_z)**2)**0.5

    return(r,r_nopbc,transVec,mic)

def calc_relVec(Ag1,com):
    pos_array = Ag1.positions
    return(pos_array - com)

residual_arrayA = np.zeros((n_frames,len(pro_A_bb)))

step = 20
chunks = np.arange(0,n_frames,step)
fr_count = 0
print(type(u.trajectory))
t1 = time_mod.perf_counter()

"""
--------------------------------------- RMSF -----------------------------------
"""



for chunk in chunks:
    u.transfer_to_memory(chunk,chunk+step)
    # print("Before: ",u.atoms.positions[100])
    protein.unwrap(reference='cog',inplace=True)
    aligner_A = align.AlignTraj(u, cry_univ,select='bynum 1-458 and name BB',in_memory=True).run()
    # print("After: ",u.atoms.positions[100])
    # print("Crystal: ",cry_univ.atoms.positions[100])
    print("Analyzing frames: ",chunk , " - ", chunk+step)
    t3 = time_mod.perf_counter()

    
    for i in range(chunk+step):
        

        comProA = pro_A.center_of_geometry(pbc=False)
        comProB = pro_B.center_of_geometry(pbc=False)

        
        # (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
        # relVec = calc_relVec(pro_B,comProB)
        # new_comProB = comProB + transVec
        # new_positions_B = relVec + new_comProB
        # if mic:
        #     print("Before: ", protein_rms.positions[469,:])

        #     pro_B_mic.translate(transVec)

        #     print("Positions after: ", protein_rms.positions[469,:])
        #     print("After chainB_mic: ", pro_B_mic.positions[0,:])

        
        residualA = distances.dist(pro_A_bb,proA_bb_cry,box=box_dims)/10.0
        # print(np.mean(residualA[2,:]))
        residual_arrayA[fr_count,:] = residualA[2,:]

        if (fr_count == chunk+step-1) or (fr_count == n_frames-1):
            fr_count+=1
            print("Break")
            break
        u.trajectory.next()
        fr_count+=1

    u.trajectory = u_copy.trajectory


print(residual_arrayA[0])
print(residual_arrayA.shape)

residual_array_sq = np.square(residual_arrayA)
avg_residuals = np.mean(residual_array_sq,axis=0)
rms_residuals = np.sqrt(avg_residuals)            #Convert to nm

rmsd_residuals = np.mean(residual_array_sq,axis=1)
rmsd_residuals = np.sqrt(rmsd_residuals)    

outfileA = outfile_pref + "_Residuals_chainA_CRYSTAL.dat"
with open(outfileA, 'w') as fl3:
    fl3.write("#Residue\tRMSF\n")
    for i in range(len(pro_A_bb.atoms)):
        fl3.write("%d\t%.3f" %(i,rms_residuals[i]))
        fl3.write("\n")

        

    fl3.write("\n")

outfileB = outfile_pref + "_RMSD_chainA_CRYSTAL.dat"
with open(outfileB,'w') as fl1:
    fl1.write("Timestep\tRMSD\n")
    for i in range(n_frames):
        fl1.write("%s\t%.3f\n" %(str(u_copy.trajectory[i].time),rmsd_residuals[i]))
        # fl1.write("\n")


"""
Old Code
"""

# aligner_B = align.AlignTraj(u_copy, cry_univ,select='bynum 470-938 and name BB',in_memory=True).run()    #The ref structure for B only contains 469 atoms of chain B. 


# rmsf_chainB_sel = u_copy.select_atoms('bynum 470-938 and name BB')
# rmsf_chainB = rms.RMSF(rmsf_chainB_sel).run()
# # rmsf_chainAtoB = rms.RMSF(rmsf_chainA_sel).run()

# outfileA = outfile_pref + "_RMSF_chainA_CRYSTAL.dat"
# with open(outfileA, 'a') as fl3:
#     fl3.write("#Residue\tRMSF\n")
#     for i in range(len(rmsf_chainA_sel.atoms)):
#         fl3.write("%d\t%.3f" %(i,rmsf_chainA.rmsf[i]))
#         fl3.write("\n")

        

#     fl3.write("\n")

# outfileB = outfile_pref + "_RMSF_chainB_CRYSTAL.dat"
# with open(outfileB, 'a') as fl3:
#     fl3.write("#Residue\tRMSF\n")
    
#     for i in range(len(rmsf_chainB_sel.atoms)):
#         fl3.write("%d\t%.3f" %(i,rmsf_chainB.rmsf[i]))
#         fl3.write("\n")

#     fl3.write("\n")