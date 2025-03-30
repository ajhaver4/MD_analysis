import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances as dist
from MDAnalysis.lib.distances import distance_array as d_array
import math
import MDAnalysis.analysis.rms as rms
import time as time_mod
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
print("Total mass: ",protein.total_mass())
# sys.exit()
#print(protein)

#Selecting protein chains
pro_A = protein.select_atoms("bynum 1-469")
pro_B = protein.select_atoms("bynum 470-938")
pro_B_mic = protein_rms.select_atoms("bynum 470-938")
atom_len = len(pro_A)
pro_A_bb = pro_A.select_atoms("name BB")
pro_B_bb = pro_B.select_atoms("name BB")

"""
Atom Group definitions for different patches
"""

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

"""
Long axis along each chain
"""
#Chain A
longA1 = protein.select_atoms("bynum 242")
#longA1 = protein.select_atoms("bynum 245")
#longA2 = protein.select_atoms("bynum 389")
longA2 = protein.select_atoms("bynum 92")

#Chain B

longB1 = protein.select_atoms("bynum 711")
#longB1 = protein.select_atoms("bynum 714")
#longB2 = protein.select_atoms("bynum 858")
longB2 = protein.select_atoms("bynum 561")

"""
Membrane binding patches
"""
memb_A = protein.select_atoms("bynum 9-11 or bynum 25-27 or bynum 32-34 or bynum 40-42 or bynum 165-167 or bynum 174-176 or bynum 181-183")
memb_B = protein.select_atoms("bynum 478-480 or bynum 494-496 or bynum 501-503 or bynum 509-511 or bynum 634-636 or bynum 643-644 or bynum 650-652")


#Reference angle b/w d1 sites:
ref_site1A = cry_univ.select_atoms("bynum 12-29") # First binding site in chain A
ref_site1B = cry_univ.select_atoms("bynum 524-556") #First binding site in chain B

"""
Defining beads for calculating pitch and yaw
"""
#Pitch
#Pick pair of beads to calculate the angle between the vector and z-axis
#Beads along the longitudnal axis

#Using only backbone beads!!
pitchA_bead1 = protein.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
pitchA_bead2 = protein.select_atoms("bynum 219") #Residue 151(pdb) 100THR (.gro)
pitchA_bead3 = protein.select_atoms("bynum 63 or bynum 70 or bynum 74")   #Residue 31GLU,35GLN,37SER (.gro)
pitchA_bead4 = protein.select_atoms("bynum 201 or bynum 212 or bynum 217")    #Residue 93ARG, 97GLU, 99ILE (.gro)

pitchB_bead1 = protein.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
pitchB_bead2 = protein.select_atoms("bynum 688") #Residue 151(pdb) 100THR (.gro)
pitchB_bead3 = protein.select_atoms("bynum 532 or bynum 539 or bynum 543")    #Residue 31GLU,35GLN,37SER (.gro)
pitchB_bead4 = protein.select_atoms("bynum 670 or bynum 681 or bynum 686")    #Residue 93ARG, 97GLU, 99ILE (.gro)

#Roll
#Chain A
base_chainA = protein.select_atoms("bynum 9")  #Residue 56 (pdb) 6ARG (.gro)
rollA_bead1 = protein.select_atoms("bynum 146")  #Residue 120 (pdb) 70ASP(.gro)
rollA_bead2 = protein.select_atoms("bynum 349")  #Residue 208(pdb) 158GLU (.gro)
rollA_bead3 = protein.select_atoms("bynum 454")  #Residue 259(pdb) 209GLU(.gro)
rollA_bead4 = protein.select_atoms("bynum 60")  #Residue 81(pdb) 30ARG (.gro)
rollA_bead5 = protein.select_atoms("bynum 9 or bynum 17")    #Residue 6ARG,9SER (.gro)
rollA_bead6 = protein.select_atoms("bynum 136 or bynum 148")    #Residue 66ASP,71LYS (.gro)

base_chainB = protein.select_atoms("bynum 478")  #Residue 56 (pdb) 6ARG (.gro)
rollB_bead1 = protein.select_atoms("bynum 615")  #Residue 120 (pdb) 70ASP(.gro)
rollB_bead2 = protein.select_atoms("bynum 818")  #Residue 208(pdb) 158GLU (.gro)
rollB_bead3 = protein.select_atoms("bynum 923")  #Residue 259(pdb) 209GLU(.gro)
rollB_bead4 = protein.select_atoms("bynum 529")  #Residue 81(pdb) 30ARG (.gro)
rollB_bead5 = protein.select_atoms("bynum 478 or bynum 486")    #Residue 6ARG,9SER (.gro)
rollB_bead6 = protein.select_atoms("bynum 605 or bynum 617")    #Residue 66ASP,71LYS (.gro)



#Printing general variables
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

def calc_angle(Ag1,Ag2,Ag3,Ag4):
    s1A = Ag1.centroid(pbc=False,compound='segments')
    s2A = Ag2.centroid(pbc=False,compound='segments')
    s1B = Ag3.centroid(pbc=False,compound='segments')
    s2B = Ag4.centroid(pbc=False,compound='segments')
    #print("Vector:",s1A)
    rA = s1A - s2A
    rB = s1B - s2B
    #print(rA)
    #print(rB)
    mod_rA = (rA[0][0]**2 + rA[0][1]**2 + rA[0][2]**2)**0.5
    mod_rB = (rB[0][0]**2 + rB[0][1]**2 + rB[0][2]**2)**0.5
    theta = math.acos((np.dot(rA.reshape(3),rB.reshape(3)))/(mod_rA*mod_rB))

    #dot_p = rA[0][0]*rB[0][0] + rA[0][1]*rB[0][1] +rA[0][2]*rB[0][2]
    #theta = math.acos((dot_p/(mod_rA*mod_rB)))

    return(theta)

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

def calc_rot(Ag1,Ag2,Ag3,Ag4,Ag5,Ag6):

    vec_A1 = Ag3.centroid(pbc=False,compound='segments') - Ag1.centroid(pbc=False,compound='segments')
    vec_A2 = Ag3.centroid(pbc=False,compound='segments') - Ag2.centroid(pbc=False,compound='segments')
    vec_B1 = Ag6.centroid(pbc=False,compound='segments') - Ag4.centroid(pbc=False,compound='segments')
    vec_B2 = Ag6.centroid(pbc=False,compound='segments') - Ag5.centroid(pbc=False,compound='segments')

    norm_A = np.cross(vec_A1,vec_A2)
    norm_B = np.cross(vec_B1,vec_B2)

    mod_rA = (norm_A[0][0]**2 + norm_A[0][1]**2 + norm_A[0][2]**2)**0.5
    mod_rB = (norm_B[0][0]**2 + norm_B[0][1]**2 + norm_B[0][2]**2)**0.5
    theta = math.acos((np.dot(norm_A.reshape(3),norm_B.reshape(3)))/(mod_rA*mod_rB))
    return(theta)

def calc_angleZ(Ag1,Ag2):
    # v1 = Ag1.centroid(pbc=False,compound='segments') - Ag2.centroid(pbc=False,compound='segments')

    v1 = Ag1.center_of_mass(pbc=False) - Ag2.center_of_mass(pbc=False)
    # print(v1)
    mod_v1 = np.sum(v1**2)**0.5
    return(math.degrees(math.acos(v1[2]/mod_v1)))

#Defining a np.array to hold the values of all variables
d1_F = np.zeros((n_frames,1))
d2_F = np.zeros((n_frames,1))
# d3_F = np.zeros((n_frames,1))
memb_dist = np.zeros((n_frames,1))
# angle_F = np.zeros((n_frames,1))
langle_F = np.zeros((n_frames,1))
# theta = np.zeros((n_frames,1))
# theta2 = np.zeros((n_frames,1))
dist_mat = np.zeros((atom_len,atom_len))
min_dist_mat = np.zeros((n_frames,1))
# memb_angle = np.zeros((n_frames,1))
rmsd_values_A = np.zeros((1,1))
rmsd_values_B = np.zeros((1,1))
#Pitch angles
pitch_A1 = np.zeros((n_frames,1))
pitch_A2 = np.zeros((n_frames,1))
pitch_A3 = np.zeros((n_frames,1))
pitch_A4 = np.zeros((n_frames,1))

pitch_B1 = np.zeros((n_frames,1))
pitch_B2 = np.zeros((n_frames,1))
pitch_B3 = np.zeros((n_frames,1))
pitch_B4 = np.zeros((n_frames,1))

#Roll Angles
roll_A1 = np.zeros((n_frames,1))
roll_A2 = np.zeros((n_frames,1))
roll_A3 = np.zeros((n_frames,1))
roll_A4 = np.zeros((n_frames,1))
roll_A5 = np.zeros((n_frames,1))

roll_B1 = np.zeros((n_frames,1))
roll_B2 = np.zeros((n_frames,1))
roll_B3 = np.zeros((n_frames,1))
roll_B4 = np.zeros((n_frames,1))
roll_B5 = np.zeros((n_frames,1))

#DeltaZ
com_z_A = np.zeros((n_frames,1))
com_z_B = np.zeros((n_frames,1))
com_z_A_membres = np.zeros((n_frames,1))
com_z_B_membres = np.zeros((n_frames,1))
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
        # print(protein.center_of_mass())
        # protein.wrap(compound='segments',center='cog',inplace=True)
        # protein_rms.wrap(compound='segments',center='cog',inplace=True)
        protein.unwrap(reference='cog',inplace=True)
        protein_rms.unwrap(reference='cog',inplace=True)
        ts = u.trajectory.ts
        # print("Timestep: ",u_copy.trajectory[fr_count].time)
        d1_F[fr_count] = calc_dist(site_1A,site_1B,False)
        d2_F[fr_count] = calc_dist(site_2A,site_2B,False)
        # d3_F[i] = calc_dist(site_3A,site_3B,False)
        memb_dist[fr_count] = calc_dist(memb_A,memb_B,False)
        # angle_F[i] = calc_angle(angle_1A,angle_2A,angle_1B,angle_2B)
        langle_F[fr_count] = math.degrees(calc_angle(longA1,longA2,longB1,longB2))
        # theta[i]=180 - math.degrees(calc_angle(hel_angle2A,hel_angle1A,hel_angle2B,hel_angle1B))
        # theta2[i] = 180 - math.degrees(calc_angle(out_hel_angle2A,out_hel_angle1A,out_hel_angle2B,out_hel_angle1B))
        # memb_angle[i] = math.degrees(calc_angle(memb_vecA,memb_A,memb_vecB,memb_B))

        # #Adjusting periodicity
        # #This is require because while the MIC is taken care while calculting individual distances,
        # #it is not considered in RMSD or dist_matrix calculations.
        # #This is problematic for those structures which are actually bound but appear unbound due to broken molecules being made whole in gmx traj
        # #Therfore the follwing code translates chain B to its minimm periodic image.
        comProA = pro_A.centroid(pbc=False)
        comProB = pro_B.centroid(pbc=False)
        (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
        relVec = calc_relVec(pro_B,comProB)
        new_comProB = comProB + transVec
        new_positions_B = relVec + new_comProB
        new_array = ts.positions
        new_array[469:938,:] = new_positions_B
        #ts._replace_positions_array(new_positions_B)
        if mic:
            print("Before: ", protein_rms.positions[469,:])
            #print("Timestep: ", ts.positions[458,:])
            # protein_rms.positions[469:938,:] = new_positions_B
            # pro_B_mic.translate(transVec)
        #     #ts._replace_positions_array(new_array)
        #     #print("Trans Vec: ",transVec)
            print("Positions after: ", protein_rms.positions[469,:])
            print("After chainB_mic: ", pro_B_mic.positions[0,:])
        #     #print("Timestep: ", ts.positions[458,:])

        dist_mat = d_array(pro_A.positions,pro_B_mic.positions,box=box_dims)
        min_dist_mat[fr_count] = np.amin(dist_mat)

        #Pitch angles
        pitch_A1[fr_count] = calc_angleZ(longA1,base_chainA)
        pitch_A2[fr_count] = calc_angleZ(longA2,base_chainA)
        pitch_A3[fr_count] = calc_angleZ(pitchA_bead1,pitchA_bead2)
        pitch_A4[fr_count] = calc_angleZ(pitchA_bead3,pitchA_bead4)

        pitch_B1[fr_count] = calc_angleZ(longB1,base_chainB)
        pitch_B2[fr_count] = calc_angleZ(longB2,base_chainB)
        pitch_B3[fr_count] = calc_angleZ(pitchB_bead1,pitchB_bead2)
        pitch_B4[fr_count] = calc_angleZ(pitchB_bead3,pitchB_bead4)

        #Roll angles
        roll_A1[fr_count] = calc_angleZ(rollA_bead1,base_chainA)
        roll_A2[fr_count] = calc_angleZ(rollA_bead2,base_chainA)
        roll_A3[fr_count] = calc_angleZ(rollA_bead3,base_chainA)
        roll_A4[fr_count] = calc_angleZ(rollA_bead3,rollA_bead4)
        roll_A5[fr_count] = calc_angleZ(rollA_bead5,rollA_bead6)

        roll_B1[fr_count] = calc_angleZ(rollB_bead1,base_chainB)
        roll_B2[fr_count] = calc_angleZ(rollB_bead2,base_chainB)
        roll_B3[fr_count] = calc_angleZ(rollB_bead3,base_chainB)
        roll_B4[fr_count] = calc_angleZ(rollB_bead3,rollB_bead4)
        roll_B5[fr_count] = calc_angleZ(rollB_bead5,rollB_bead6)

        #DeltaZ
        com_z_A[fr_count]=comProA[2]
        com_z_B[fr_count] = comProB[2]

        com_z_A_membres[fr_count] = memb_A.center_of_geometry(pbc=False)[2] 
        com_z_B_membres[fr_count] = memb_B.center_of_geometry(pbc=False)[2] 


        if (fr_count == chunk+step-1) or (fr_count == n_frames-1):
            fr_count+=1
            print("Break")

            break
        u.trajectory.next()
        fr_count+=1
    # #Resetting the trajectory for RMSD analysis
    """
    --------------------------------------- RMSD -----------------------------------
    """

    lsp_rmsd= rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 1-469', groupselections=["name BB and bynum 1-469", "name BB and (bynum 524-556 or bynum 770-788)"])
    lsp_rmsd.run()
    size_arr = lsp_rmsd.rmsd.T[3].shape[0]
    print("Size: ",size_arr)

    rmsd_values_A = np.concatenate((rmsd_values_A,np.reshape(lsp_rmsd.rmsd.T[3],(size_arr,1))),axis=0)
    rmsd_values_B = np.concatenate((rmsd_values_B,np.reshape(lsp_rmsd.rmsd.T[4],(size_arr,1))),axis=0)

    u.trajectory = u_copy.trajectory
    t4 = time_mod.perf_counter()
    print("Time taken for reading chunk %.4f" %(t4-t3))

u.trajectory.rewind()
t2 = time_mod.perf_counter()
print("Time taken for complete analysis: %.4f" %(t2-t1))


#Delete first rows
rmsd_values_A = np.delete(rmsd_values_A,0,0)
rmsd_values_B = np.delete(rmsd_values_B,0,0)

"""
#--------------------------------------------------------------------------------
#--------------------------- PLOTTING/WRITING TRAJ -------------------------------------------
#--------------------------------------------------------------------------------
"""
time_str = 'Simualtion time: ' + str(time)
output_file = "../Final_paramter_data.txt"
with open(output_file,'a') as fl1:
    fl1.write('##')
    fl1.write('\t')
    fl1.write(time_str)
    fl1.write('\n')
    fl1.write("#Frame No.")
    fl1.write("\t")
    fl1.write("Timestep")
    fl1.write("\t")
    fl1.write("d1")
    fl1.write("\t")
    fl1.write("d2")
    fl1.write("\t")
    # fl1.write("dist")
    # fl1.write("\t")
    # fl1.write("ang")
    # fl1.write("\t")
    fl1.write("Ax-ang")
    fl1.write("\t")
    # fl1.write("Dim_ang-1")
    # fl1.write("\t")
    # fl1.write("Dim_ang-2")
    # fl1.write("\t")
    fl1.write("RMSD-A")
    fl1.write("\t")
    fl1.write("RMSD-B")
    fl1.write("\t")
    fl1.write("Min dist")
    fl1.write("\t")
    fl1.write("Memb dist")
    # fl1.write("\t")
    # fl1.write("Rot_angle")
    fl1.write("\n")

    for i in range(n_frames):
        fl1.write(str(i))
        fl1.write("\t")
        fl1.write(str(u_copy.trajectory[i].time))
        fl1.write("\t")
        fl1.write(str(d1_F[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(d2_F[i][0]/10.0))
        fl1.write("\t")
        # fl1.write(str(d3_F[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(angle_F[i][0]))
        # fl1.write("\t")
        fl1.write(str(langle_F[i][0]))
        fl1.write("\t")
        # fl1.write(str(theta[i][0]))
        # fl1.write("\t")
        # fl1.write(str(theta2[i][0]))
        # fl1.write("\t")
        fl1.write(str(rmsd_values_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_B[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(min_dist_mat[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(memb_dist[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(memb_angle[i][0]))
        fl1.write("\n")


with open("../Delta_Z_measurements.txt",'a') as fl1:
    # fl1.write("#Frame\tTimestep\tcomA\tcomB\tdelta_z_A\tdelta_z_B\tz_A_hi\tz_B_hi\tz_A_low\tz_B_low\tz_A_med\tz_B_med\tz_A_membres\tz_B_membres\n")
    fl1.write("#Frame\tTimestep\tcomA\tcomB\tz_A_membres\tz_B_membres\n")
    for i in range(n_frames):
        fl1.write(str(i))
        fl1.write("\t")
        fl1.write(str(u_copy.trajectory[i].time))
        fl1.write("\t")
        fl1.write(str(com_z_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_B[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_A_membres[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_B_membres[i][0]/10.0))
        fl1.write("\n")

with open("../Pitch_Roll_measurements.txt", 'a') as fl2:
    fl2.write("#Frame\tTimestep\tPitchA1\tPitchA2\tPitchA3\tPitchA4\tPitchB1\tPitchB2\tPitchB3\tPitchB4\tRollA1\tRollA2\tRollA3\tRollA4\tRollA5\tRollB1\tRollB2\tRollB3\tRollB4\tRollB5\n")
    for i in range(n_frames):
        fl2.write(str(i))
        fl2.write("\t")
        fl2.write(str(u_copy.trajectory[i].time))
        fl2.write("\t")
        fl2.write(str(pitch_A1[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A2[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A3[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A4[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_B1[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_B2[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_B3[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_B4[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_A1[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_A2[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_A3[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_A4[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_A5[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_B1[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_B2[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_B3[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_B4[i][0]))
        fl2.write("\t")
        fl2.write(str(roll_B5[i][0]))
        fl2.write("\n")





