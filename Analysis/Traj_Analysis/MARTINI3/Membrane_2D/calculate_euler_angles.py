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
from scipy.spatial.transform import Rotation
from MDAnalysis.lib.transformations import *
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


#Selecting protein chains
pro_A = protein.select_atoms("bynum 1-469")
pro_B = protein.select_atoms("bynum 470-938")
pro_B_mic = protein_rms.select_atoms("bynum 470-938")
atom_len = len(pro_A)
pro_A_bb = pro_A.select_atoms("name BB")
pro_B_bb = pro_B.select_atoms("name BB")


proA_bb_cry = cry_univ.select_atoms("bynum 1-469 and name BB")
proB_bb_cry = cry_univ.select_atoms("bynum 470-938 and name BB")

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
longA1 = protein.select_atoms("bynum 242")     #Residues 159(.pdb) 109LYS(.gro)
#longA1 = protein.select_atoms("bynum 245")
#longA2 = protein.select_atoms("bynum 389")
longA2 = protein.select_atoms("bynum 92")       #Residues 95 (.pdb) 45ASP (.gro)

#Chain B

longB1 = protein.select_atoms("bynum 711")
#longB1 = protein.select_atoms("bynum 714")
#longB2 = protein.select_atoms("bynum 858")
longB2 = protein.select_atoms("bynum 561")

"""
Membrane binding patches
"""
#6ARG , 13 LYS, 16 LYS, 20ARG, 76ARG, 80LYS, 83ARG  -> (.gro)
#56ARG , 63LYS, 66LYS, 70ARG, 126ARG, 130LYS, 133ARG -> (.pdb)
#Chain A
memb_A = protein.select_atoms("bynum 9-11 or bynum 25-27 or bynum 32-34 or bynum 40-42 or bynum 165-167 or bynum 174-176 or bynum 181-183")
#Chain B
memb_B = protein.select_atoms("bynum 478-480 or bynum 494-496 or bynum 501-503 or bynum 509-511 or bynum 634-636 or bynum 643-644 or bynum 650-652")



"""
Defining beads for calculating pitch and yaw
"""
#Pitch
#Pick pair of beads to calculate the angle between the vector and z-axis
#Beads along the longitudnal axis

#Using only backbone beads!!
pitchA_bead1 = protein.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
pitchA_bead2 = protein.select_atoms("bynum 219") #Residue 151(pdb) 100THR (.gro)

pitchA_bead1_cry = cry_univ.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
pitchA_bead2_cry = cry_univ.select_atoms("bynum 219") #Residue 151(pdb) 100THR (.gro)

pitchB_bead1 = protein.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
pitchB_bead2 = protein.select_atoms("bynum 688") #Residue 151(pdb) 100THR (.gro)

pitchB_bead1_cry = cry_univ.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
pitchB_bead2_cry = cry_univ.select_atoms("bynum 688") #Residue 151(pdb) 100THR (.gro)

#Roll
#Chain A
rollA_bead3 = protein.select_atoms("bynum 454")  #Residue 259(pdb) 209GLU(.gro)
rollA_bead4 = protein.select_atoms("bynum 60")  #Residue 81(pdb) 30ARG (.gro)

rollB_bead3 = protein.select_atoms("bynum 923")  #Residue 259(pdb) 209GLU(.gro)
rollB_bead4 = protein.select_atoms("bynum 529")  #Residue 81(pdb) 30ARG (.gro)


#Beads for defining long axis (axis of rotation for Roll), perp axis (axis of rotation for yaw),  short axis (axis of rotation for pitch)

#For the long-axis we will use the same definition as the beads used for pitch 
#For defining perp axis, one bead is the same as the one used for Pitch, the other is a slightly different bead (close to the second roll bead). Could have used the same, but the new bead seems more perpendicular.
bead_perpA1 = protein.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
bead_perpA2 = protein.select_atoms("bynum 447") #Residue 255(pdb) 205GLU (.gro)

bead_perpA1_cry = cry_univ.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
bead_perpA2_cry = cry_univ.select_atoms("bynum 447") #Residue 255(pdb) 205GLU (.gro)

bead_perpB1 = protein.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
bead_perpB2 = protein.select_atoms("bynum 916") #Residue 255(pdb) 205GLU (.gro)

bead_perpB1_cry = cry_univ.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
bead_perpB2_cry = cry_univ.select_atoms("bynum 916") #Residue 255(pdb) 205GLU (.gro)

#For the short axis, again keeping one bead the same as the one used for pitch and roll (Kind of like the origin or rotation). 
bead_shortA1 = protein.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
bead_shortA2 = protein.select_atoms("bynum 12") #Residues 57 (pdb) 7LYS (.gro)

bead_shortA1_cry = cry_univ.select_atoms("bynum 67") #Residue 84(pdb) 34LYS(.gro)
bead_shortA2_cry = cry_univ.select_atoms("bynum 12") #Residues 57 (pdb) 7LYS (.gro)

bead_shortB1 = protein.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
bead_shortB2 = protein.select_atoms("bynum 481") #Residues 57 (pdb) 7LYS (.gro)

bead_shortB1_cry = cry_univ.select_atoms("bynum 536") #Residue 84(pdb) 34LYS(.gro)
bead_shortB2_cry = cry_univ.select_atoms("bynum 481") #Residues 57 (pdb) 7LYS (.gro)

#Reference angle b/w d1 sites:
ref_site1A = cry_univ.select_atoms("bynum 12-29") # First binding site in chain A
ref_site1B = cry_univ.select_atoms("bynum 524-556") #First binding site in chain B
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

def calc_angleZ(Ag1,Ag2):
    # v1 = Ag1.centroid(pbc=False,compound='segments') - Ag2.centroid(pbc=False,compound='segments')

    v1 = Ag1.center_of_mass(pbc=False) - Ag2.center_of_mass(pbc=False)
    # print(v1)
    mod_v1 = np.sum(v1**2)**0.5
    return(math.degrees(math.acos(v1[2]/mod_v1)))

def calc_Yaw(Ag1,Ag2):
    # v1 = Ag1.centroid(pbc=False,compound='segments') - Ag2.centroid(pbc=False,compound='segments')

    v1 = Ag1.center_of_mass(pbc=False) - Ag2.center_of_mass(pbc=False)
    
    angle = math.atan2(v1[1],v1[0])

    if angle < 0:
        angle = 2*math.pi + angle
    return(math.degrees(angle))

def calc_eulers(l_axis,l_perp):

    l_axis = l_axis/np.linalg.norm(l_axis)
    l_perp = l_perp/np.linalg.norm(l_perp)

    z1 = np.array([0,0,1])
    pitch_axis = np.cross(z1,l_axis)
    pitch_axis = pitch_axis/np.linalg.norm(pitch_axis)
    yaw_axis = np.cross(l_axis,pitch_axis)

    
    #Angles
    psi_elem = np.arctan2(l_axis[1],l_axis[0])
    theta_elem = np.arcsin(-l_axis[2])
    phi_elem = np.arctan2(np.dot(l_perp,yaw_axis),np.dot(l_perp,pitch_axis))


    # init_angle = math.atan2(np.inner(init_vec,refvec2),np.inner(init_vec,refvec1)) 
    # curr_angle = math.atan2(np.inner(curr_vec,refvec2),np.inner(curr_vec,refvec1))

    # rel_roll = curr_angle - init_angle
    # if phi_elem < 0:
    #     phi_elem = 2*math.pi + phi_elem

    # if psi_elem < 0:
    #     psi_elem = 2*math.pi + psi_elem

    return([psi_elem,theta_elem,phi_elem])

def calc_eulerMg(l_axis,l_perp):
    l_axis = l_axis/np.linalg.norm(l_axis)
    l_perp = l_perp/np.linalg.norm(l_perp)

    pitch_axis = np.cross(l_perp,l_axis)
    pitch_axis = pitch_axis/np.linalg.norm(pitch_axis)
    yaw_axis = np.cross(l_axis,pitch_axis)

    
    #Angles
    psi_elem = np.arctan2(l_axis[1],l_axis[0])
    theta_elem = np.arcsin(-l_axis[2])
    phi_elem = np.arctan2(pitch_axis[2],yaw_axis[2])


    return([psi_elem,theta_elem,phi_elem])



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
rmsd_values_BtoA = np.zeros((1,1))
rmsd_values_BtoB = np.zeros((1,1))

pitchA_angle = np.zeros((n_frames,1))
rollA_angle = np.zeros((n_frames,1))
rollB_angle = np.zeros((n_frames,1))
pitchB_angle = np.zeros((n_frames,1))

yawA_angle = np.zeros((n_frames,1))
yawB_angle = np.zeros((n_frames,1))


#Angles calculated from rotation matrices
rotMat_pitchA = np.zeros((n_frames,1))
rotMat_pitchB = np.zeros((n_frames,1))
rotMat_rollA = np.zeros((n_frames,1))
rotMat_rollB = np.zeros((n_frames,1))
rotMat_yawA = np.zeros((n_frames,1))
rotMat_yawB = np.zeros((n_frames,1))


pitchA_angleMg = np.zeros((n_frames,1))
rollA_angleMg = np.zeros((n_frames,1))
rollB_angleMg = np.zeros((n_frames,1))
pitchB_angleMg = np.zeros((n_frames,1))
yawA_angleMg = np.zeros((n_frames,1))
yawB_angleMg = np.zeros((n_frames,1))





#Defining 2 vectors for the reference frame which is perpendicular to the Pitch vector



init_rlongA = pitchA_bead2.center_of_mass() - pitchA_bead1.center_of_mass()    #Doesn't matter if pbs True or not since it is just the first frame used a a reference vector
init_rlongB = pitchB_bead2.center_of_mass() - pitchB_bead1.center_of_mass()

init_rshortA = bead_perpA2.center_of_mass() - bead_perpA1.center_of_mass()
init_rshortB = bead_perpB2.center_of_mass() - bead_perpB1.center_of_mass()


init_rollAngleA = calc_eulers(init_rlongA,init_rshortA)[2]
init_rollAngleB = calc_eulers(init_rlongB,init_rshortB)[2]

print("Initial Roll Angle A: ", math.degrees(init_rollAngleA))
print("Initial Roll Angle B: ", math.degrees(init_rollAngleB))

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
        protein.wrap(compound='segments',center='cog',inplace=True)
        protein_rms.wrap(compound='segments',center='cog',inplace=True)
        # protein.unwrap(reference='cog',inplace=True)
        # protein_rms.unwrap(reference='cog',inplace=True)
        ts = u.trajectory.ts
        d1_F[fr_count] = calc_dist(site_1A,site_1B,False)
        d2_F[fr_count] = calc_dist(site_2A,site_2B,False)
        
        memb_dist[fr_count] = calc_dist(memb_A,memb_B,False)
        
        langle_F[fr_count] = math.degrees(calc_angle(longA1,longA2,longB1,longB2))
        

        # #Adjusting periodicity
        # #This is require because while the MIC is taken care while calculting individual distances,
        # #it is not considered in RMSD or dist_matrix calculations.
        # #This is problematic for those structures which are actually bound but appear unbound due to broken molecules being made whole in gmx traj
        # #Therfore the follwing code translates chain B to its minimm periodic image.
        comProA = pro_A.centroid(pbc=False,compound='segments')
        comProB = pro_B.centroid(pbc=False,compound='segments')
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

        #Calculate Pitch angle
        rlong_vecA = pitchA_bead2.center_of_mass() - pitchA_bead1.center_of_mass()
        rlong_vecB = pitchB_bead2.center_of_mass() - pitchB_bead1.center_of_mass()

        rshort_vecA = bead_perpA2.center_of_mass() - bead_perpA1.center_of_mass()  #rollA_bead3.center_of_mass() - rollA_bead4.center_of_mass()
        rshort_vecB = bead_perpB2.center_of_mass() - bead_perpB1.center_of_mass()  #rollB_bead3.center_of_mass() - rollB_bead4.center_of_mass()

        curr_anglesA = calc_eulers(rlong_vecA,rshort_vecA)     #Returns angles in radians
        curr_anglesB = calc_eulers(rlong_vecB,rshort_vecB)

        pitchA_angle[fr_count] = curr_anglesA[1]
        pitchB_angle[fr_count] = curr_anglesB[1]

        rollA_angle[fr_count] = curr_anglesA[2] #- init_rollAngleA
        rollB_angle[fr_count] = curr_anglesB[2] #- init_rollAngleB

        yawA_angle[fr_count] = curr_anglesA[0]
        yawB_angle[fr_count] = curr_anglesB[0]

        #Mg method
        curr_anglesA_Mg = calc_eulerMg(rlong_vecA,rshort_vecA)     #Returns angles in radians
        curr_anglesB_Mg = calc_eulerMg(rlong_vecB,rshort_vecB)

        pitchA_angleMg[fr_count] = curr_anglesA_Mg[1]
        pitchB_angleMg[fr_count] = curr_anglesB_Mg[1]
        rollA_angleMg[fr_count] = curr_anglesA_Mg[2]
        rollB_angleMg[fr_count] = curr_anglesB_Mg[2]
        yawA_angleMg[fr_count] = curr_anglesA_Mg[0]
        yawB_angleMg[fr_count] = curr_anglesB_Mg[0]

        #ALigning trajectory for calculating Euler angles

        #Two types of alignment. First calculating using backbone atoms. Then using atoms that define only principle axis. 
        print("Starting alignment for chain A")
        print("COM of chain B: ", pro_B.centroid(pbc=False,compound='segments'))
        # rotMatrixA, rmsdA = align.rotation_matrix(pro_A_bb.positions,proA_bb_cry.positions)
        # r = Rotation.from_matrix(rotMatrixA)
        # euler_anglesA = r.as_euler('xyz', degrees=True)

        # rotMat_yawA[fr_count] = euler_anglesA[0]
        # rotMat_pitchA[fr_count] = euler_anglesA[1]
        # rotMat_rollA[fr_count] = euler_anglesA[2]

        align_posbeads = np.concatenate((pitchA_bead1.positions,pitchA_bead2.positions,bead_perpA1.positions,bead_perpA2.positions,bead_shortA1.positions,bead_shortA2.positions))
        ref_posbeads = np.concatenate((pitchA_bead1_cry.positions,pitchA_bead2_cry.positions,bead_perpA1_cry.positions,bead_perpA2_cry.positions,bead_shortA1_cry.positions,bead_shortA2_cry.positions))
        rotMatrixA, rmsdA = align.rotation_matrix(ref_posbeads,align_posbeads)
        r = Rotation.from_matrix(rotMatrixA)
        euler_anglesA = r.as_euler('xyz', degrees=True)

        rotMat_yawA[fr_count] = euler_anglesA[0]
        rotMat_pitchA[fr_count] = euler_anglesA[1]
        rotMat_rollA[fr_count] = euler_anglesA[2]
        

        al_out = align.alignto(u,cry_univ,select='bynum 1-469 and name BB')

        print("Initializing rotation for chain B: ")
        print("COM of chain B: ", pro_B.centroid(pbc=False,compound='segments'))
        rotMatrixB,rmsdB = align.rotation_matrix(pro_B_bb.positions,proB_bb_cry.positions)
        rB = Rotation.from_matrix(rotMatrixB)
        euler_anglesB = rB.as_euler('xyz', degrees=True)

        rotMat_yawB[fr_count] = euler_anglesB[0]
        rotMat_pitchB[fr_count] = euler_anglesB[1]
        rotMat_rollB[fr_count] = euler_anglesB[2]

        
        print("COM of B: ", pro_B.centroid(pbc=False,compound='segments'))

        if (fr_count == chunk+step-1) or (fr_count == n_frames-1):
            fr_count+=1
            print("Break")

            break
        u.trajectory.next()
        fr_count+=1
    
    


    u.trajectory = u_copy.trajectory
    t4 = time_mod.perf_counter()
    print("Time taken for reading chunk %.4f" %(t4-t3))

u.trajectory.rewind()
t2 = time_mod.perf_counter()
print("Time taken for complete analysis: %.4f" %(t2-t1))



"""
#--------------------------------------------------------------------------------
#--------------------------- PLOTTING/WRITING TRAJ -------------------------------------------
#--------------------------------------------------------------------------------
"""
time_str = 'Simualtion time: ' + str(time)
output_file = "EulerAngles_data.txt"
with open(output_file,'a') as fl2:
    fl2.write("#Frame\tTimestep\tPitchA\tPitchB\tRollA\tRollB\tYawA\tYawB\trotMatPitchA\trotMatPitchB\trotMatRollA\trotMatRollB\trotMatYawA\trotMatYawB\n")
    for i in range(n_frames):
        fl2.write(str(i))
        fl2.write("\t")
        fl2.write(str(u_copy.trajectory[i].time))
        fl2.write("\t")
        fl2.write(str(pitchA_angle[i][0]))
        fl2.write("\t")
        fl2.write(str(pitchB_angle[i][0]))
        fl2.write("\t")
        fl2.write(str(rollA_angle[i][0]))
        fl2.write("\t")
        fl2.write(str(rollB_angle[i][0]))
        fl2.write("\t")
        fl2.write(str(yawA_angle[i][0]))
        fl2.write("\t")
        fl2.write(str(yawB_angle[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_pitchA[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_pitchB[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_rollA[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_rollB[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_yawA[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_yawB[i][0]))
        fl2.write("\n")


output_file = "EulerAngles_dataMg.txt"
with open(output_file,'a') as fl2:
    fl2.write("#Frame\tTimestep\tPitchA\tPitchB\tRollA\tRollB\tYawA\tYawB\trotMatPitchA\trotMatPitchB\trotMatRollA\trotMatRollB\trotMatYawA\trotMatYawB\n")
    for i in range(n_frames):
        fl2.write(str(i))
        fl2.write("\t")
        fl2.write(str(u_copy.trajectory[i].time))
        fl2.write("\t")
        fl2.write(str(pitchA_angleMg[i][0]))
        fl2.write("\t")
        fl2.write(str(pitchB_angleMg[i][0]))
        fl2.write("\t")
        fl2.write(str(rollA_angleMg[i][0]))
        fl2.write("\t")
        fl2.write(str(rollB_angleMg[i][0]))
        fl2.write("\t")
        fl2.write(str(yawA_angleMg[i][0]))
        fl2.write("\t")
        fl2.write(str(yawB_angleMg[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_pitchA[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_pitchB[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_rollA[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_rollB[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_yawA[i][0]))
        fl2.write("\t")
        fl2.write(str(rotMat_yawB[i][0]))
        fl2.write("\n")



