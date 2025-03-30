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
time = sys.argv[4]

#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()
#Selecting only protein beads
protein = u.select_atoms("byres name BB")
protein_rms = u.select_atoms("byres name BB")
#print(protein)

#Selecting protein chains
pro_A = u.select_atoms("bynum 1-458")
pro_B = u.select_atoms("bynum 459-916")
pro_B_mic = protein_rms.select_atoms("bynum 459-916")
atom_len = len(pro_A)
pro_A_bb = pro_A.select_atoms("name BB")
pro_B_bb = pro_B.select_atoms("name BB")

#Defining Atom Groups for dist_CV
#These groups correspond to the binding sites for dist_CV in metadynamics
#Chains A & B
#Sites 1 and 2 in each chain
site_1A = protein.select_atoms("bynum 12-29") # First binding site in chain A
site_1B = protein.select_atoms("bynum 513-544") #First binding site in chain B
site_2A = protein.select_atoms("bynum 375-386") #Second binding site in chain A
site_2B = protein.select_atoms("bynum 754-772") #Second binding site in chain B

#Defining Atom Groups for angle_CV
site_3A = protein.select_atoms("bynum 337 or bynum 338 or bynum 344 or bynum 345 or bynum 346 or bynum 347 or bynum 353 or bynum 354 or bynum 355 or bynum 356")
site_3B = protein.select_atoms("bynum 795 or bynum 796 or bynum 802 or bynum 803 or bynum 804 or bynum 805 or bynum 811 or bynum 812 or bynum 813 or bynum 814")

angle_1A = protein.select_atoms("bynum 19-20")
angle_2A = protein.select_atoms("bynum 444-445")
angle_1B = protein.select_atoms("bynum 477-478")
angle_2B = protein.select_atoms("bynum 902-903")

#Define angles along the axis of molecules
langle_1A = protein.select_atoms("bynum 244")
langle_2A = protein.select_atoms("bynum 395")

langle_1B = protein.select_atoms("bynum 702")
langle_2B = protein.select_atoms("bynum 853")

#Definning AtomGroup for Angle of Dimerisation
#theta1
#End atoms (two residues) of Helix 2 in each chain
hel_angle1A = protein.select_atoms("name BB and bynum 6-8")
hel_angle2A = protein.select_atoms("name BB and bynum 233-239")
hel_angle1B = protein.select_atoms("name BB and bynum 464-466")
hel_angle2B = protein.select_atoms("name BB and bynum 691-697")

#theta2
#Considering only the outer arms of helix2
out_hel_angle1A = protein.select_atoms("name BB and bynum 189-192")
out_hel_angle2A = protein.select_atoms("name BB and bynum 233-239")
out_hel_angle1B = protein.select_atoms("name BB and bynum 647-650")
out_hel_angle2B = protein.select_atoms("name BB and bynum 691-697")

#Membrane binding sites
memb_A = protein.select_atoms("bynum 9-11 or bynum 25-27 or bynum 32-34 or bynum 40-42 or bynum 161-163 or bynum 170-172 or bynum 177-179")
memb_B = protein.select_atoms("bynum 467-469 or bynum 483-485 or bynum 490-492 or bynum 498-500 or bynum 619-621 or bynum 628-630 or bynum 635-637")

memb_vecA = protein.select_atoms("bynum 305 or bynum 307 or bynum 316 or bynum 323")
memb_vecB = protein.select_atoms("bynum 763 or bynum 820 or bynum 823 or bynum 833 or bynum 765 or bynum 779")
#Reference angle b/w d1 sites:
ref_site1A = cry_univ.select_atoms("bynum 12-29") # First binding site in chain A
ref_site1B = cry_univ.select_atoms("bynum 513-544") #First binding site in chain B
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
#Defining a np.array to hold the values of all variables
d1_F = np.zeros((n_frames,1))
d2_F = np.zeros((n_frames,1))
d3_F = np.zeros((n_frames,1))
memb_dist = np.zeros((n_frames,1))
angle_F = np.zeros((n_frames,1))
langle_F = np.zeros((n_frames,1))
theta = np.zeros((n_frames,1))
theta2 = np.zeros((n_frames,1))
dist_mat = np.zeros((atom_len,atom_len))
min_dist_mat = np.zeros((n_frames,1))
memb_angle = np.zeros((n_frames,1))
rmsd_values_A = np.zeros((1,1))
rmsd_values_B = np.zeros((1,1))
#Looping over each frame to calculate d1 and d2 for each frame and the store them in the above arrays
#u.transfer_to_memory()
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
        d3_F[fr_count] = calc_dist(site_3A,site_3B,False)
        memb_dist[fr_count] = calc_dist(memb_A,memb_B,False)
        angle_F[fr_count] = calc_angle(angle_1A,angle_2A,angle_1B,angle_2B)
        langle_F[fr_count] = calc_angle(langle_1A,langle_2A,langle_1B,langle_2B)
        theta[fr_count]=180 - math.degrees(calc_angle(hel_angle2A,hel_angle1A,hel_angle2B,hel_angle1B))
        theta2[fr_count] = 180 - math.degrees(calc_angle(out_hel_angle2A,out_hel_angle1A,out_hel_angle2B,out_hel_angle1B))
        memb_angle[fr_count] = math.degrees(calc_angle(memb_vecA,memb_A,memb_vecB,memb_B))


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
        new_array[458:916,:] = new_positions_B
        #ts._replace_positions_array(new_positions_B)
        if mic:
            print("Before: ", protein_rms.positions[458,:])
            #print("Timestep: ", ts.positions[458,:])
            # protein_rms.positions[469:938,:] = new_positions_B
            # pro_B_mic.translate(transVec)
        #     #ts._replace_positions_array(new_array)
        #     #print("Trans Vec: ",transVec)
            print("Positions after: ", protein_rms.positions[458,:])
            print("After chainB_mic: ", pro_B_mic.positions[0,:])
        #     #print("Timestep: ", ts.positions[458,:])

        dist_mat = d_array(pro_A.positions,pro_B_mic.positions,box=box_dims)
        min_dist_mat[fr_count] = np.amin(dist_mat)

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

    lsp_rmsd= rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 1-458', groupselections=["name BB and bynum 1-458", "name BB and (bynum 513-544 or bynum 754-772)"])
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
output_file = "Final_data.txt"
with open(output_file,'a') as fl1:
    fl1.write('##')
    fl1.write('\t')
    fl1.write(time_str)
    fl1.write('\n')
    fl1.write("#Frame No. ")
    fl1.write("\t")
    fl1.write("Timestep")
    fl1.write("\t")
    fl1.write("d1")
    fl1.write("\t")
    fl1.write("d2")
    fl1.write("\t")
    fl1.write("dist")
    fl1.write("\t")
    fl1.write("ang")
    fl1.write("\t")
    fl1.write("Ax-ang")
    fl1.write("\t")
    fl1.write("Dim_ang-1")
    fl1.write("\t")
    fl1.write("Dim_ang-2")
    fl1.write("\t")
    fl1.write("RMSD-A")
    fl1.write("\t")
    fl1.write("RMSD-B")
    fl1.write("\t")
    fl1.write("Min dist")
    fl1.write("\t")
    fl1.write("Memb dist")
    fl1.write("\t")
    fl1.write("Rot_angle")
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
        fl1.write(str(d3_F[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(angle_F[i][0]))
        fl1.write("\t")
        fl1.write(str(langle_F[i][0]))
        fl1.write("\t")
        fl1.write(str(theta[i][0]))
        fl1.write("\t")
        fl1.write(str(theta2[i][0]))
        fl1.write("\t")
        fl1.write(str(rmsd_values_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_B[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(min_dist_mat[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(memb_dist[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(memb_angle[i][0]))
        fl1.write("\n")
