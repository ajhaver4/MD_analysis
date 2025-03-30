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
"""
Inputs - struct_file traj_file num_of_clusters crystal_struct_file
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-------------------------------- INITIALIZATION --------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
"""
#Reading in all the inputs
struct_file = sys.argv[1]
traj_file = sys.argv[2]
cry_univ = mda.Universe(sys.argv[3],sys.argv[4])
time = sys.argv[5]

#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()

#Selecting only protein beads
protein = u.select_atoms("byres name BB")
protein_rms = u.select_atoms("byres name BB")

chainA = protein.select_atoms("bynum 1-469")
chainB = protein.select_atoms("bynum 470-938")
chainB_mic = protein_rms.select_atoms("bynum 470-938")

#Selecting protein chains
pro_A = protein.select_atoms("bynum 1-469")
pro_B = protein.select_atoms("bynum 470-938")



n_frames = u.trajectory.n_frames
ts = u.trajectory.ts
print("Timestep:",ts.dt)
global box_dims
box_dims = ts.dimensions
print("Box dimension: ", box_dims)
print("Number of frames: ", n_frames)

"""
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
------------------------- Defining Methods ------------------------------------
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
"""

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

def calc_dist(Ag1,Ag2,pbc):
    global box_dims
    s1A = Ag1.centroid(pbc=pbc,compound='segments')
    s1B = Ag2.centroid(pbc=pbc,compound='segments')
    #mic_val = np.zeros((3))
    #del_nopbc = np.zeros((3))
    #del_pbc = np.zeros((3))
    #print(len(Ag1.residues))
    #print(len(Ag2.residues))
    #print(s1A.shape)
    #print(s1B.shape)
    #print("First coords:", s1A)
    #print("Second coords:", s1B)

    del_x = s1A[0][0] - s1B[0][0]
    del_y = s1A[0][1] - s1B[0][1]
    del_z = s1A[0][2] - s1B[0][2]

    del_x = del_x - (box_dims[0])*round(del_x/box_dims[0],0)
    del_y = del_y - (box_dims[1])*round(del_y/box_dims[1],0)
    del_z = del_z - (box_dims[2])*round(del_z/box_dims[2],0)


    r = ((del_x)**2 + (del_y)**2 + (del_z)**2)**0.5
    #dist = distance.euclidean(s1A,s1B)
    return(r)
#Defining a np.array to hold the values of all variables
"""
METAD CV patches
"""
#d1
#Chains A & B
#Sites 1 and 2 in each chain
site_1A = protein.select_atoms("bynum 12-29") # First binding site in chain A
site_1B = protein.select_atoms("bynum 524-556") #First binding site in chain B
site_2A = protein.select_atoms("bynum 383-394") #Second binding site in chain A
site_2B = protein.select_atoms("bynum 770-788") #Second binding site in chain B

comZ_chainA = np.zeros((n_frames,1))
comZ_chainB  = np.zeros((n_frames,1))

rmsd_values_A = np.zeros((1,1))
rmsd_values_BtoA = np.zeros((1,1))
rmsd_values_BtoB = np.zeros((1,1))

bound_frames = []
bound_ts = []

unbound_frames = []
unbound_ts = []

#Looping over each frame to calculate d1 and d2 for each frame and the store them in the above arrays
step = 20
chunks = np.arange(0,n_frames,step)
fr_count = 0
print(type(u.trajectory))
t1 = time_mod.perf_counter()
for chunk in chunks:
    u.transfer_to_memory(chunk,chunk+step)
    # print(type(u.trajectory))
    # u.trajectory = u_copy.trajectory
    # print(type(u.trajectory))
    # sys.exit()
    print("Analyzing frames: ",chunk , " - ", chunk+step)
    t3 = time_mod.perf_counter()
    for i in range(chunk+step):

        protein.unwrap(reference='cog',inplace=True)
        ts = u.trajectory.ts
        curr_time = u_copy.trajectory[fr_count].time
        comProA = pro_A.centroid(pbc=False,compound='segments')
        comProB = pro_B.centroid(pbc=False,compound='segments')

        # #Check which protein has crossed the boundary
        # out_protein = check_out_of_bounds(comProA,comProB)
        # if out_protein=='B':
        #     (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
        #     relVec = calc_relVec(chainB,comProB)
        #     new_comProB = comProB + transVec
        #     new_positions_B = relVec + new_comProB
        # if out_protein=='A':
        #     (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProB,3),np.reshape(comProA,3))
        #     relVec = calc_relVec(chainA,comProA)
        #     new_comProA = comProA + transVec
        #     new_positions_A = relVec + new_comProA
        (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
        relVec = calc_relVec(chainB,comProB)
        new_comProB = comProB + transVec
        new_positions_B = relVec + new_comProB
        if mic:
            print("Before: ", protein_rms.positions[470,:])
            #print("Timestep: ", ts.positions[458,:])
            # protein_rms.positions[469:938,:] = new_positions_B
            chainB_mic.translate(transVec)
        #     #ts._replace_positions_array(new_array)
        #     #print("Trans Vec: ",transVec)
            print("Positions after: ", protein_rms.positions[470,:])
            print("After chainB_mic: ", chainB_mic.positions[0,:])
        
        
        # print("Timestep: ",u_copy.trajectory[fr_count].time)
        curr_d1 = calc_dist(site_1A,site_1B,False)/10.0
        curr_d2 = calc_dist(site_2A,site_2B,False)/10.0
        #Storing comZ values
        comZ_chainA[fr_count] = comProA[0][2]
        comZ_chainB[fr_count] = comProB[0][2]

        if curr_d1 < 13.0 and curr_d2 < 13.0:
            bound_frames.append(fr_count)
            bound_ts.append(curr_time)
        else:
            unbound_frames.append(fr_count)
            unbound_ts.append(curr_time)


        if (fr_count == chunk+step-1) or (fr_count == n_frames-1):
            fr_count+=1
            print("Break")

            break
        u.trajectory.next()
        fr_count+=1
    

    """
    --------------------------------------- RMSD -----------------------------------
    """

    lsp_rmsd= rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 1-469', groupselections=["name BB and bynum 1-469", "name BB and (bynum 524-556 or bynum 770-788)"])
    lsp_rmsd.run()
    size_arr = lsp_rmsd.rmsd.T[3].shape[0]
    print("Size: ",size_arr)

    rmsd_wrtB = rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 470-938')
    rmsd_wrtB.run()

    rmsd_values_A = np.concatenate((rmsd_values_A,np.reshape(lsp_rmsd.rmsd.T[3],(size_arr,1))),axis=0)    #The first two colums are frame and timesteps, the 3rd is the rmsd for the superimposed group.
    rmsd_values_BtoA = np.concatenate((rmsd_values_BtoA,np.reshape(lsp_rmsd.rmsd.T[4],(size_arr,1))),axis=0)
    rmsd_values_BtoB = np.concatenate((rmsd_values_BtoB,np.reshape(rmsd_wrtB.rmsd.T[2],(size_arr,1))),axis=0)
    
    u.trajectory = u_copy.trajectory
    t4 = time_mod.perf_counter()
    print("Time taken for reading chunk %.4f" %(t4-t3))
    

    


#Delete first rows
rmsd_values_A = np.delete(rmsd_values_A,0,0)
rmsd_values_BtoA = np.delete(rmsd_values_BtoA,0,0)
rmsd_values_BtoB = np.delete(rmsd_values_BtoB,0,0)

u.trajectory.rewind()
t2 = time_mod.perf_counter()
print("Time taken for complete analysis: %.4f" %(t2-t1))


bd_lbl = "traj_bound_" + str(time) + ".xtc"
unbd_lbl = "traj_unbound_" + str(time) + ".xtc"

protein.write(bd_lbl,frames=bound_frames)
protein.write(unbd_lbl,frames=unbound_frames)

print("Finished writing bound and unbound trajectories")

"""
#--------------------------------------------------------------------------------
#--------------------------- PLOTTING/WRITING TRAJ -------------------------------------------
#--------------------------------------------------------------------------------
"""
time_str = 'Simualtion time: ' + str(time)
output_file = "COM_data.txt"
with open(output_file,'a') as fl1:
    fl1.write('##')
    fl1.write('\t')
    fl1.write(time_str)
    fl1.write('\n')
    fl1.write("#Frame No.")
    fl1.write("\t")
    fl1.write("Timestep")
    fl1.write("\t")
    fl1.write("COMz_chainA")
    fl1.write("\t")
    fl1.write("COMz_chainB")
    fl1.write("\t")
    fl1.write("RMSD-A")
    fl1.write("\t")
    fl1.write("RMSD-BtoA")
    fl1.write("\t")
    fl1.write("RMSD-BtoB")
    fl1.write("\n")

    for i in range(n_frames):
        fl1.write(str(i))
        fl1.write("\t")
        fl1.write(str(u_copy.trajectory[i].time))
        fl1.write("\t")
        fl1.write(str(comZ_chainA[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(comZ_chainB[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_BtoA[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_BtoB[i][0]/10.0))
        fl1.write("\n")



