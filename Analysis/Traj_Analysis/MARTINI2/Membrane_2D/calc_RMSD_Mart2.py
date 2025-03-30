import sys
#import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import math
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import distance_array as d_array
from MDAnalysis.analysis.contacts import contact_matrix
from MDAnalysis.core.groups import AtomGroup
import time as time_mod
import MDAnalysis.analysis.rms as rms
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt

struct_file = sys.argv[1]
traj_file = sys.argv[2]
cry_univ = mda.Universe(sys.argv[3],sys.argv[4])


#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()
#Protein selcetions
protein=u.select_atoms("byres name BB")     #All protein beads (1094)
protein_rms = u.select_atoms("byres name BB")

chainA = protein.select_atoms("bynum 1-458")
chainB = protein.select_atoms("bynum 459-916")
chainB_mic = protein_rms.select_atoms("bynum 459-916")


rmsd_values_A = np.zeros((1,1))
rmsd_values_BtoA = np.zeros((1,1))
rmsd_values_BtoB = np.zeros((1,1))

n_frames = u.trajectory.n_frames
ts = u.trajectory.ts
global box_dims
box_dims = ts.dimensions

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

step = 200
chunks = np.arange(0,n_frames,step)
fr_count = 0
print(type(u.trajectory))
t1 = time_mod.perf_counter()
for chunk in chunks:
    u.transfer_to_memory(chunk,chunk+step)
    print("Analyzing frames: ",chunk , " - ", chunk+step)
    t3 = time_mod.perf_counter()

    
    for i in range(chunk+step):
        protein.unwrap(reference='cog',inplace=True)

        comProA = chainA.center_of_geometry(pbc=False)
        comProB = chainB.center_of_geometry(pbc=False)
        
        (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
        relVec = calc_relVec(chainB,comProB)
        new_comProB = comProB + transVec
        new_positions_B = relVec + new_comProB
        if mic:
            print("Before: ", protein_rms.positions[470,:])

            chainB_mic.translate(transVec)

            print("Positions after: ", protein_rms.positions[470,:])
            print("After chainB_mic: ", chainB_mic.positions[0,:])


        if (fr_count == chunk+step-1) or (fr_count == n_frames-1):
            fr_count+=1
            print("Break")
            break
        u.trajectory.next()
        fr_count+=1

    """
    --------------------------------------- RMSD -----------------------------------
    """

    lsp_rmsd= rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 1-458', groupselections=["name BB and bynum 1-458", "name BB and (bynum 513-544 or bynum 754-772)"])
    lsp_rmsd.run()
    size_arr = lsp_rmsd.rmsd.T[3].shape[0]
    print("Size: ",size_arr)

    rmsd_wrtB = rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 459-916')
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
t2 = time_mod.perf_counter()
print("Time taken for complete analysis: %.4f" %(t2-t1))

outfile= "RMSD_wrt_" + str(sys.argv[4]) + "_.dat"
with open(outfile,"a") as fl:
    fl.write("#Frame\tTime\tRMSD-A\tRMSD-BtoA\tRMSD-BtoB\n")
    for i in range(n_frames):
        fl.write("%d\t%s\t%.3f\t%.3f\t%.3f\n" %(i,str(u_copy.trajectory[i].time),rmsd_values_A[i][0]/10.0,rmsd_values_BtoA[i][0]/10.0,rmsd_values_BtoB[i][0]/10.0))
        
       