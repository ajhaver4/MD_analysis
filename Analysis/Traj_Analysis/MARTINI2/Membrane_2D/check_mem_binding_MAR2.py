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

struct_file = sys.argv[1]
traj_file = sys.argv[2]
cry_univ = mda.Universe(sys.argv[3],sys.argv[4])
time = sys.argv[5]

#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()
#Protein selcetions
protein=u.select_atoms("byres name BB")     #All protein beads (1094)
protein_rms = u.select_atoms("byres name BB")

chainA = protein.select_atoms("bynum 52425-52882")
chainB = protein.select_atoms("bynum 52883-53340")
chainB_mic = protein_rms.select_atoms("bynum 52883-53340")

side_chains = protein.select_atoms("(byres name BB) and not name BB")  #Only side chain beads (590)
backbone = protein.select_atoms("name BB")    #Only backbone beads(504)
print(chainA)
#sys.exit()
#Select membrane atoms
memb = u.select_atoms("byres name PO4")     #All lipid beads (54728); Not a useful selection
#General parameters
po4 = u.select_atoms("name PO4")        #Selecting phospho head groups for all lipids (4414)
print(po4)
print(len(po4))
# sys.exit()
popc_lipid = u.select_atoms("resname POPC and name PO4") #Selecting phospho headgroups for diff lipids (3534)
pops_lipid = u.select_atoms("resname POPS and name PO4") #(440)
pop6_lipid = u.select_atoms("resname POP6 and name PO4") #(440)

com_po4 = po4.center_of_geometry()
upp_leaf_ind = []
low_leaf_ind = []
for i in range(len(po4)):
    z = po4[i].position[2]
    if z > com_po4[2]:
        upp_leaf_ind.append(i)
    elif z < com_po4[2]:
        low_leaf_ind.append(i)


upp_leaf = po4[upp_leaf_ind]
low_leaf = po4[low_leaf_ind]

# print(upp_leaf.positions[2454][2])
# print(low_leaf[0].position[2])
# sys.exit()
"""
Atom Group definitions for different patches
"""

"""
METAD CV patches
"""
#d1
#Chains A & B
#Sites 1 and 2 in each chain
#In this system, residue numbering of protein starts from 60797 (Atom 1 on chainA) and ends at  161872 (Atom 2 on chainB)
#Atom no. in this sytem is obtained by Atom no. in solution system - 9 + 60797  -> Chain A
#Atom no. in this sytem is obtained by Atom no. in solution system - 18 + 60797  -> Chain B

site_1A = protein.select_atoms("bynum 52436-52453") # First binding site in chain A         #Residues 57LYS-> 64THR  (pdb numbering)
site_1B = protein.select_atoms("bynum 52937-52968") #First binding site in chain B              #Residues 78GLU-92ASP
site_2A = protein.select_atoms("bynum 52799-52810") #Second binding site in chain A             #Residues 224LEU-229ASP
site_2B = protein.select_atoms("bynum 53178-53196") #Second binding site in chain B             #Residues 189ASN-196LYS

"""
Long axis along each chain
"""
#Chain A
#longA1 = protein.select_atoms("name BB and resid 180 and segid seg_0_Protein_A")
longA1 = protein.select_atoms("bynum 52661")
#longA1 = protein.select_atoms("bynum 245")
#longA2 = protein.select_atoms("name BB and resid 252 and segid A")
longA2 = protein.select_atoms("bynum 52819")
#longA2 = protein.select_atoms("name BB and resid 114 and segid seg_0_Protein_A")

#Chain B

#longB1 = protein.select_atoms("name BB and resid 432 and segid seg_1_Protein_B")
longB1 = protein.select_atoms("bynum 53119")
#longB1 = protein.select_atoms("bynum 714")
#longB2 = protein.select_atoms("name BB and resid 252 and segid B")
#longB2 = protein.select_atoms("name BB and resid 366 and segid seg_1_Protein_B")
longB2 = protein.select_atoms("bynum 53277")

"""
Membrane binding patches
"""
memb_A = protein.select_atoms("bynum 52433-52435 or bynum 52449-52451 or bynum 52456-52458 or bynum 52464-52466 or bynum 52594-52596 or bynum 52601-52603")
memb_B = protein.select_atoms("bynum 52891-52893 or bynum 52907-52909 or bynum 52914-52916 or bynum 52922-52924 or bynum 53052-53054 or bynum 53059-53061")

"""
Dimerisation Angle
"""
##Definning AtomGroup for Angle of Dimerisation
#theta1
#End atoms (two residues) of Helix 2 in each chain
hel_angle1A = protein.select_atoms("name BB and bynum 52430-52432")
hel_angle2A = protein.select_atoms("name BB and bynum 52657-52663")
hel_angle1B = protein.select_atoms("name BB and bynum 52888-52890")
hel_angle2B = protein.select_atoms("name BB and bynum 53115-53121")

#theta2
#Considering only the outer arms of helix2
out_hel_angle1A = protein.select_atoms("name BB and bynum 52613-52616")
out_hel_angle2A = protein.select_atoms("name BB and bynum 52657-52663")
out_hel_angle1B = protein.select_atoms("name BB and bynum 53071-53074")
out_hel_angle2B = protein.select_atoms("name BB and bynum 53115-53121")


"""
End-to-End distance patches
"""
end_chainA1 = protein.select_atoms("bynum 52819")
end_chainA2 = protein.select_atoms("bynum 52661")

end_chainB1 = protein.select_atoms("bynum 53277")
end_chainB2 = protein.select_atoms("bynum 53119")


"""
Restraining Wall patches
"""
# wall_A = protein.select_atoms("name BB and resid 5021-5024 and segid seg_620_Protein_P")
#wall_A = protein.select_atoms("name BB and resid 112-115")
# wall_B = protein.select_atoms("name BB and resid 5262-5265 and segid seg_621_Protein_Q")
#wall_B = protein.select_atoms("name BB and resid 364-367")


"""
Defining beads to calculate COM of each chain
"""
com_A_hi_flx = protein.select_atoms("bynum 52657 or bynum 52661 or bynum 52664 or bynum 52666 or bynum 52668 or bynum 52672")
com_B_hi_flx = protein.select_atoms("bynum 53115 or bynum 53117 or bynum 53119 or bynum 53122 or bynum 53125 or bynum 53126")
#
com_A_low_flx = protein.select_atoms("bynum 52527 or bynum 52533 or bynum 52535 or bynum 52538 or bynum 52540 or bynum 52542")
com_B_low_flx = protein.select_atoms("bynum 52984 or bynum 52988 or bynum 53991 or bynum 53993 or bynum 53996 or bynum 53998")
#
com_A_med_flx = protein.select_atoms("bynum 52657 or bynum 52533 or bynum 52664 or bynum 52538 or bynum 52668 or bynum 52542")
com_B_med_flx = protein.select_atoms("bynum 53115 or bynum 52988 or bynum 53119 or bynum 53993 or bynum 53125 or bynum 53998")


"""
Defining beads for calculating pitch and yaw
"""
#Pitch
#Pick pair of beads to calculate the angle between the vector and z-axis
#Beads along the longitudnal axis

pitchA_bead1 = protein.select_atoms("bynum 52491") #Residue 84(pdb) 34LYS(.gro)
pitchA_bead2 = protein.select_atoms("bynum 52639") #Residue 151(pdb) 100THR (.gro)
pitchA_bead3 = protein.select_atoms("bynum 52487 or bynum 52494 or bynum 52498")   #Residue 31GLU,35GLN,37SER (.gro)
pitchA_bead4 = protein.select_atoms("bynum 52621 or bynum 52632 or bynum 52637")    #Residue 93ARG, 97GLU, 99ILE (.gro)

pitchB_bead1 = protein.select_atoms("bynum 52949") #Residue 84(pdb) 34LYS(.gro)
pitchB_bead2 = protein.select_atoms("bynum 53097") #Residue 151(pdb) 100THR (.gro)
pitchB_bead3 = protein.select_atoms("bynum 52945 or bynum 52952 or bynum 52956")    #Residue 31GLU,35GLN,37SER (.gro)
pitchB_bead4 = protein.select_atoms("bynum 53079 or bynum 53090 or bynum 53095")    #Residue 93ARG, 97GLU, 99ILE (.gro)

#Roll
#Chain A
base_chainA = protein.select_atoms("bynum 52433")  #Residue 56 (pdb) 6ARG (.gro)
rollA_bead1 = protein.select_atoms("bynum 52568")  #Residue 120 (pdb) 70ASP(.gro)
rollA_bead2 = protein.select_atoms("bynum 52766")  #Residue 208(pdb) 158GLU (.gro)
rollA_bead3 = protein.select_atoms("bynum 52868")  #Residue 259(pdb) 209GLU(.gro)
rollA_bead4 = protein.select_atoms("bynum 52484")  #Residue 81(pdb) 30ARG (.gro)
rollA_bead5 = protein.select_atoms("bynum 52433 or bynum 52441")    #Residue 6ARG,9SER (.gro)
rollA_bead6 = protein.select_atoms("bynum 52558 or bynum 52570")    #Residue 66ASP,71LYS (.gro)

base_chainB = protein.select_atoms("bynum 52891")  #Residue 56 (pdb) 6ARG (.gro)
rollB_bead1 = protein.select_atoms("bynum 53026")  #Residue 120 (pdb) 70ASP(.gro)
rollB_bead2 = protein.select_atoms("bynum 53224")  #Residue 208(pdb) 158GLU (.gro)
rollB_bead3 = protein.select_atoms("bynum 53326")  #Residue 259(pdb) 209GLU(.gro)
rollB_bead4 = protein.select_atoms("bynum 52942")  #Residue 81(pdb) 30ARG (.gro)
rollB_bead5 = protein.select_atoms("bynum 52891 or bynum 52899")    #Residue 6ARG,9SER (.gro)
rollB_bead6 = protein.select_atoms("bynum 53016 or bynum 53028")    #Residue 66ASP,71LYS (.gro)



# print(len(upp_leaf))
# print(len(low_leaf))

# print(len(memb))
# print(len(popc_lipid))
# print(len(pops_lipid))
# print(len(pop6_lipid))
# print(len(po4))
# print(len(protein))
# print(len(side_chains))
# print(len(backbone))
#sys.exit()
ts = u.trajectory.ts
global box_dims
box_dims = ts.dimensions
#Find COGs of each group
# cog_protein = protein.center_of_geometry(compound='segments',unwrap=True)
#The unwrap argument can be tricky. If using for delta_Z calculations double check if pbc conditions are met even when unwrap is False
#For other calculations unwrap can be False. On second thoughts since protein is stuck on membrane the z_coordinate will face any pbc issues
try:
    cog_protein = protein.center_of_geometry(compound='segments',unwrap=True)
except:
    cog_protein = protein.center_of_geometry(compound='segments',unwrap=False)

cog_upp = upp_leaf.center_of_geometry()
cog_low = low_leaf.center_of_geometry()
cog_po4 = po4.center_of_geometry()

print(cog_protein)
print(cog_upp)
print(cog_low)
# sys.exit()
print(chainA.center_of_geometry())
print(chainB.center_of_geometry())


# sys.exit()
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
    # s1A = Ag1.center_of_geometry(pbc=False)
    # s2A = Ag2.center_of_geometry(pbc=False)
    # s1B = Ag3.center_of_geometry(pbc=False)
    # s2B = Ag4.center_of_geometry(pbc=False)

    s1A = Ag1.centroid(pbc=False,compound='segments')
    s2A = Ag2.centroid(pbc=False,compound='segments')
    s1B = Ag3.centroid(pbc=False,compound='segments')
    s2B = Ag4.centroid(pbc=False,compound='segments')
    #print("Vector:",s1B)
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

def check_out_of_bounds(comA,comB):
    out_pro = ''
    for i in range(len(comA)):
        if comA[i] < 0 or comA[i] > box_dims[i]:
            out_pro='A'
    for j in range(len(comB)):
        if comB[j] < 0 or comB[j] > box_dims[j]:
            out_pro='B'

    return(out_pro)

def calc_angleZ(Ag1,Ag2):
    # v1 = Ag1.centroid(pbc=False,compound='segments') - Ag2.centroid(pbc=False,compound='segments')

    v1 = Ag1.center_of_mass(pbc=False) - Ag2.center_of_mass(pbc=False)
    # print(v1)
    mod_v1 = np.sum(v1**2)**0.5
    return(math.degrees(math.acos(v1[2]/mod_v1)))
# sys.exit()
##Check for closest atoms for each side chain bead
#Find distances between protien sidechains and phospho headgroups
#depending upon whether the protein is closer to upper leaflet or lower leaflet

class Protein_Data():
    def __init__(self,group_len,n_frames):
        self.total_count = np.zeros((group_len,1))
        self.bound_frames = 0
        self.membraneBoundbool = False
        self.boundIntervals = []
        self.inTime = 0.0
        self.outTime = 0.0
        self.totalFrames = 0.0

n_frames = u.trajectory.n_frames
min_dist_lipid = np.zeros((n_frames,1))
min_dist_A_lipid = np.zeros((n_frames,1))
min_dist_B_lipid = np.zeros((n_frames,1))

memb_min_distA_lipid = np.zeros((n_frames,1))
memb_min_distB_lipid = np.zeros((n_frames,1))
com_protein = np.zeros((n_frames,1))
min_residues = []
mic_array = []
contact_mat = np.zeros((len(side_chains),len(upp_leaf)))

n_frames_bound = 0

#Defining distance variables
d1_F = np.zeros((n_frames,1))
d2_F = np.zeros((n_frames,1))
l_angle = np.zeros((n_frames,1))
dim_angle = np.zeros((n_frames,1))

min_dist_mat_pro = np.zeros((n_frames,1))

# ins_angle = np.zeros((n_frames,1))
memb_dist = np.zeros((n_frames,1))
com_z_A = np.zeros((n_frames,1))
com_z_B = np.zeros((n_frames,1))

com_z_A_hi_flx = np.zeros((n_frames,1))
com_z_B_hi_flx = np.zeros((n_frames,1))
com_z_A_low_flx = np.zeros((n_frames,1))
com_z_B_low_flx = np.zeros((n_frames,1))
com_z_A_med_flx = np.zeros((n_frames,1))
com_z_B_med_flx = np.zeros((n_frames,1))
com_z_A_membres = np.zeros((n_frames,1))
com_z_B_membres = np.zeros((n_frames,1))

Z_angle_A = np.zeros((n_frames,1))
Z_angle_B = np.zeros((n_frames,1))

delta_z_A = np.zeros((n_frames,1))
delta_z_B = np.zeros((n_frames,1))

dz_xy_A_upp = np.zeros((n_frames,1))
dz_xy_B_upp = np.zeros((n_frames,1))

dz_xy_A_low = np.zeros((n_frames,1))
dz_xy_B_low = np.zeros((n_frames,1))

lipid_count = np.zeros((n_frames,1))
leaf_flag = ""

wall_z_A = np.zeros((n_frames,1))
wall_z_B = np.zeros((n_frames,1))

num_sc = int(len(side_chains)/2)

rmsd_values_A = np.zeros((1,1))
rmsd_values_BtoA = np.zeros((1,1))
rmsd_values_BtoB = np.zeros((1,1))

end_to_end_A = np.zeros((n_frames,1))
end_to_end_B = np.zeros((n_frames,1))

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

#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
pick_path = "./total_contacts_data.pickle"
if not os.path.exists(pick_path):
    curr_data = Protein_Data(len(side_chains),n_frames)
else:
    with open(pick_path,'rb') as pick_handle:
        curr_data = pickle.load(pick_handle)

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
        po4.unwrap(reference='cog',inplace=True)

        comProA = chainA.center_of_geometry(pbc=False)
        comProB = chainB.center_of_geometry(pbc=False)

        new_comProA = comProA
        new_comProB = comProB




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
        #protein.unwrap(compound='segments',reference='cog')
        #cog_protein = protein.center_of_geometry(compound='segments')
        chainA_comZ = comProA[2]
        chainB_comZ = comProB[2]
        #cog_upp = upp_leaf.center_of_geometry()
        #cog_low = low_leaf.center_of_geometry()
        cog_po4 = po4.center_of_geometry()

        #Calculating distance between COG of each protein chain with all lipid headgroups of each layer. Using only XY coordinates
        #The concatenation is done since the array shape has to be of (3,) or (N,3) i.e. have 3 coordinates for the d_array function input. 
        # So the Z coordinate is set to zero for the d_array calculation since we only want the min XY distance 
        xyproj_comA = np.concatenate([comProA[0:2],np.zeros((1,))])
        xyproj_comB = np.concatenate([comProB[0:2],np.zeros((1,))])

        xyproj_uppleaf = np.concatenate([upp_leaf.positions[:,0:2],np.zeros((len(upp_leaf),1))],axis=1)
        xyproj_lowleaf = np.concatenate([low_leaf.positions[:,0:2],np.zeros((len(low_leaf),1))],axis=1)

        memb_xydist_comA_upp = d_array(xyproj_comA,xyproj_uppleaf)
        memb_xydist_comA_low = d_array(xyproj_comA,xyproj_lowleaf)

        memb_xydist_comB_upp = d_array(xyproj_comB,xyproj_uppleaf)
        memb_xydist_comB_low = d_array(xyproj_comB,xyproj_lowleaf)

        #Calculating all distances between proteinA and bilayer.
        if chainA_comZ - cog_po4[2] > 0:
            # dist_array = d_array(side_chains.positions,upp_leaf.positions)
            memb_dist_array_A = d_array(memb_A.positions,upp_leaf.positions)
            leaf_flagA = 'up'
        elif chainA_comZ - cog_po4[2] < 0:
            # dist_array = d_array(side_chains.positions,low_leaf.positions)
            memb_dist_array_A = d_array(memb_A.positions,low_leaf.positions)
            leaf_flagA = 'low'

        #Calculating all distances between proteinB and bilayer.
        if chainB_comZ - cog_po4[2] > 0:
            memb_dist_array_B = d_array(memb_B.positions,upp_leaf.positions)
            leaf_flagB = 'up'
        elif chainB_comZ - cog_po4[2] < 0:
            memb_dist_array_B = d_array(memb_B.positions,low_leaf.positions)
            leaf_flagB = 'low'

        #Calculating all distances between protein and bilayer
        #Initially this was done for one particular leafleet, but since protein chains sometimes unbind and bind to other leaflet, 
        #it is better to calculate distances for both leaflets and then choose the minimum distance
        dist_array = d_array(side_chains.positions,po4.positions)

        #Getting the minimum residues for corresponding distances
        min_dist_res = np.unravel_index(np.argmin(dist_array,axis=None),dist_array.shape)
        minIndx_memb_A = np.unravel_index(np.argmin(memb_dist_array_A,axis=None),memb_dist_array_A.shape)
        minIndx_memb_B = np.unravel_index(np.argmin(memb_dist_array_B,axis=None),memb_dist_array_B.shape)
        min_des_res_chainA = np.unravel_index(np.argmin(dist_array[0:num_sc,:],axis=None),dist_array[0:num_sc,:].shape)
        min_des_res_chainB = np.unravel_index(np.argmin(dist_array[num_sc:,:],axis=None),dist_array[num_sc:,:].shape)

        #Storing the distances
        min_dist_lipid[fr_count] = dist_array[min_dist_res[0],min_dist_res[1]]/10.0
        min_dist_A_lipid[fr_count] = dist_array[min_des_res_chainA[0],min_des_res_chainA[1]]/10.0
        min_dist_B_lipid[fr_count] = dist_array[min_des_res_chainB[0]+num_sc,min_des_res_chainB[1]]/10.0
        memb_min_distA_lipid[fr_count] = memb_dist_array_A[minIndx_memb_A[0],minIndx_memb_A[1]]/10.0
        memb_min_distB_lipid[fr_count] = memb_dist_array_B[minIndx_memb_B[0],minIndx_memb_B[1]]/10.0

        # print("-------------------   Frame no. %d -------------------" %(i))
        # print("Minimum Residues: ")
        # print("Dimer: ", side_chains[min_dist_res[0]].ix,side_chains[min_dist_res[0]].resname,side_chains[min_dist_res[0]].segid)
        # print("ChainA: ",side_chains[min_des_res_chainA[0]].ix,side_chains[min_des_res_chainA[0]].resname,side_chains[min_des_res_chainA[0]].segid)
        # print("ChainB: ",side_chains[min_des_res_chainB[0]+num_sc].ix,side_chains[min_des_res_chainB[0]+num_sc].resname,side_chains[min_des_res_chainB[0]+num_sc].segid)

        #Getting the residues which are closest in xy distance to the protein COG
        xyIndx_uppA = np.unravel_index(np.argmin(memb_xydist_comA_upp,axis=None),memb_xydist_comA_upp.shape)
        xyIndx_lowA = np.unravel_index(np.argmin(memb_xydist_comA_low,axis=None),memb_xydist_comA_low.shape)
        xyIndx_uppB = np.unravel_index(np.argmin(memb_xydist_comB_upp,axis=None),memb_xydist_comB_upp.shape)
        xyIndx_lowB = np.unravel_index(np.argmin(memb_xydist_comB_low,axis=None),memb_xydist_comB_low.shape)

        #Storing the distances
        dz_xy_A_upp[fr_count] = comProA[2] - upp_leaf.positions[xyIndx_uppA[0],2]
        dz_xy_A_low[fr_count] = comProA[2] - low_leaf.positions[xyIndx_lowA[0],2]
        dz_xy_B_upp[fr_count] = comProB[2] - upp_leaf.positions[xyIndx_uppB[0],2]
        dz_xy_B_low[fr_count] = comProB[2] - low_leaf.positions[xyIndx_lowB[0],2]

        #Pick 10 closest headgroups to calculate delta_z distance
        #Taking minmum along the entire protein chains
        lipid_distances = dist_array.min(axis=0)/10.0
        print("Lipid distances: ",lipid_distances)
        #print(len(lipid_distances))
        lipid_groups_z = []
        for l in range(len(lipid_distances)):
            if lipid_distances[l] < 1.2:
                #Select this group to be in calculation for delta_z
                #Want to be sure that the lipid headgroups using for delta_z calculations
                #are closest to protein
                # if leaf_flag == 'up':
                #     lipid_groups_z.append(upp_leaf[l].position[2])
                # elif leaf_flag =='low':
                #     lipid_groups_z.append(low_leaf[l].position[2])
                lipid_groups_z.append(po4[l].position[2])

        lipid_count[fr_count] = len(lipid_groups_z)
        COM_lipid_z = np.mean(lipid_groups_z)

        print("Lipid Groups: ",len(lipid_groups_z))
        print("COM_lipid_z: ",COM_lipid_z)
        #Pick 10 closest headgroups from membrane residues to cal delz_membrane
        #Taking minmum along the entire protein chains
        lipid_distances_z_A = memb_dist_array_A.min(axis=0)/10.0
        lipid_distances_z_B = memb_dist_array_B.min(axis=0)/10.0
        #print(len(lipid_distances))
        lipid_groups_z_A = []
        lipid_groups_z_B = []
        for l in range(len(lipid_distances_z_A)):
            if lipid_distances_z_A[l] < 1.2:
                #Select this group to be in calculation for delta_z
                #Want to be sure that the lipid headgroups using for delta_z calculations
                #are closest to protein
                if leaf_flagA == 'up':
                    lipid_groups_z_A.append(upp_leaf[l].position[2])
                elif leaf_flagA =='low':
                    lipid_groups_z_A.append(low_leaf[l].position[2])

            if lipid_distances_z_B[l] < 1.2:
                #Select this group to be in calculation for delta_z
                #Want to be sure that the lipid headgroups using for delta_z calculations
                #are closest to protein
                if leaf_flagB == 'up':
                    lipid_groups_z_B.append(upp_leaf[l].position[2])
                elif leaf_flagB =='low':
                    lipid_groups_z_B.append(low_leaf[l].position[2])

        # print("Lipid groups: ",len(lipid_groups_z_A),len(lipid_groups_z_B))
        COM_lipid_z_A = np.mean(lipid_groups_z_A)
        COM_lipid_z_B = np.mean(lipid_groups_z_B)
        #print(COM_lipid_z)
        # print(COM_lipid_z_A)
        # print(COM_lipid_z_B)



        # print("Minimum distance: ",min_dist[i])
        # print("Minimum distance 2: ",np.amin(dist_array))
        # if leaf_flag=='up':
            # print("Leaflet: ",leaf_flag)
            # print("Lipid residues: ", upp_leaf[min_dist_res[1]].ix,upp_leaf[min_dist_res[1]].resname)
        # elif leaf_flag=='low':
        #     print("Leaflet: ",leaf_flag)
        #     print("Lipid residues: ", low_leaf[min_dist_res[1]].ix,low_leaf[min_dist_res[1]].resname)
        if min_dist_lipid[fr_count]<0.8:
            if not curr_data.membraneBoundbool:
                curr_data.inTime = u.trajectory[i].time
                curr_data.membraneBoundbool = True
                print("Membrane Binding Event")
                print("Time: ",u.trajectory[i].time)
                print("Residue :",side_chains[min_dist_res[0]].ix,side_chains[min_dist_res[0]].resname)
            curr_data.bound_frames+=1
        elif min_dist_lipid[fr_count]>=0.8 and curr_data.membraneBoundbool:
            curr_data.outTime = u.trajectory[i].time
            curr_data.membraneBoundbool=False
            curr_data.boundIntervals.append([curr_data.inTime,curr_data.outTime])


        com_z_A[fr_count]=comProA[2]
        com_z_B[fr_count] = comProB[2]

        #Calculting delta_Z values for different bead combinations
        delta_z_A[fr_count] = abs(comProA[2] - COM_lipid_z)
        delta_z_B[fr_count] = abs(comProB[2] - COM_lipid_z)

        com_z_A_hi_flx[fr_count] = abs(com_A_hi_flx.center_of_geometry(pbc=False)[2] - COM_lipid_z)
        com_z_B_hi_flx[fr_count] = abs(com_B_hi_flx.center_of_geometry(pbc=False)[2] - COM_lipid_z)
        com_z_A_low_flx[fr_count] = abs(com_A_low_flx.center_of_geometry(pbc=False)[2] - COM_lipid_z)
        com_z_B_low_flx[fr_count] = abs(com_B_low_flx.center_of_geometry(pbc=False)[2] - COM_lipid_z)
        com_z_A_med_flx[fr_count] = abs(com_A_med_flx.center_of_geometry(pbc=False)[2] - COM_lipid_z)
        com_z_B_med_flx[fr_count] = abs(com_B_med_flx.center_of_geometry(pbc=False)[2] - COM_lipid_z)

        com_z_A_membres[fr_count] = abs(memb_A.center_of_geometry(pbc=False)[2] - COM_lipid_z_A)
        com_z_B_membres[fr_count] = abs(memb_B.center_of_geometry(pbc=False)[2] - COM_lipid_z_B)


        min_residues.append(min_dist_res)
        # com_protein[fr_count] = protein_z
        mic_array.append(mic)

        #Contact MAtrix
        # contact_mat = contact_matrix(dist_array,8,out=contact_mat)

        contact_mat = contact_matrix(dist_array,8)   #Removing a pre-allocated matrix since at some timesteps protein is binding to diff leaflets that changes the size of  dist_array and throws an error in this line while concatenating.        
        protein_contacts = list(set(np.nonzero(contact_mat)[0]))

        for j in range(len(protein_contacts)):
            curr_data.total_count[protein_contacts[j]]+=1
        # if mic:
        #     print("Frame: ",i)
        #     print(comProA)
        #     print(comProB)
        #     print(new_comProB)
        #     print(com_protein[i-1])
        #     sys.exit()

        d1_F[fr_count] = calc_dist(site_1A,site_1B,False)/10.0
        d2_F[fr_count] = calc_dist(site_2A,site_2B,False)/10.0
        l_angle[fr_count] = math.degrees(calc_angle(longA1,longA2,longB1,longB2))
        dim_angle[fr_count] = 180 - math.degrees(calc_angle(hel_angle2A,hel_angle1A,hel_angle2B,hel_angle1B))
        memb_dist[fr_count] = calc_dist(memb_A,memb_B,False)/10.0

        dist_mat_pro = d_array(chainA.positions,chainB_mic.positions,box=box_dims)
        min_dist_mat_pro[fr_count] = np.amin(dist_mat_pro)/10.0

        end_to_end_A[fr_count] = calc_dist(end_chainA1,end_chainA2,False)/10.0
        end_to_end_B[fr_count] = calc_dist(end_chainB1,end_chainB2,False)/10.0

        Z_angle_A[fr_count] = calc_angleZ(longA1,longA2)
        Z_angle_B[fr_count] = calc_angleZ(longB1,longB2)

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


        # wall_z_A[fr_count] = wall_A.centroid(pbc=False,compound='segments')[0][2]
        # wall_z_B[fr_count] = wall_B.centroid(pbc=False,compound='segments')[0][2]

        if (fr_count == chunk+step-1) or (fr_count == n_frames-1):
            fr_count+=1
            print("Break")
            break
        u.trajectory.next()
        fr_count+=1

    """
    --------------------------------------- RMSD -----------------------------------
    """

    lsp_rmsd= rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 52425-52882', groupselections=["name BB and bynum 52425-52882", "name BB and (bynum 52937-52968 or bynum 53178-53196)"])
    lsp_rmsd.run()
    size_arr = lsp_rmsd.rmsd.T[3].shape[0]
    print("Size: ",size_arr)

    rmsd_wrtB = rms.RMSD(protein_rms,cry_univ,select='name BB and bynum 52883-53340')
    rmsd_wrtB.run()

    rmsd_values_A = np.concatenate((rmsd_values_A,np.reshape(lsp_rmsd.rmsd.T[3],(size_arr,1))),axis=0)
    rmsd_values_BtoA = np.concatenate((rmsd_values_BtoA,np.reshape(lsp_rmsd.rmsd.T[4],(size_arr,1))),axis=0)
    rmsd_values_BtoB = np.concatenate((rmsd_values_BtoB,np.reshape(rmsd_wrtB.rmsd.T[2],(size_arr,1))),axis=0)



    u.trajectory = u_copy.trajectory
    t4 = time_mod.perf_counter()
    print("Time taken for reading chunk %.4f" %(t4-t3))

t2 = time_mod.perf_counter()
print("Time taken for complete analysis: %.4f" %(t2-t1))

#Delete first rows
rmsd_values_A = np.delete(rmsd_values_A,0,0)
rmsd_values_BtoA = np.delete(rmsd_values_BtoA,0,0)
rmsd_values_BtoB = np.delete(rmsd_values_BtoB,0,0)



print("Maximum Z distances: ")
print("Chain A: ",max(wall_z_A))
print("Chain B: ",max(wall_z_B))


with open("./Final_parameter_data.txt",'a') as fl1:
    fl1.write("#Frame\tTimestep\td1\td2\tlangle\tdim-angle\tMinProDist\tmin_dist\tmin_dist_A\tmin_dist_B\tmemb_minA\tmemb_minB\tRMSD-A\tRMSD-B\tRMSD-BtoA\tcomA\tcomB\tmembdist\tzangle_A\tzangle_B\tdelta_z_A\tdelta_z_B\tlipid_count\tEnd2ENdA\tEnd2EndB\n")
    for i in range(n_frames):
        fl1.write(str(i))
        fl1.write("\t")
        fl1.write(str(u_copy.trajectory[i].time))
        fl1.write("\t")
        fl1.write(str(d1_F[i][0]))
        fl1.write("\t")
        fl1.write(str(d2_F[i][0]))
        fl1.write("\t")
        fl1.write(str(l_angle[i][0]))
        fl1.write("\t")
        fl1.write(str(dim_angle[i][0]))
        fl1.write("\t")
        fl1.write(str(min_dist_mat_pro[i][0]))
        fl1.write("\t")
        fl1.write(str(min_dist_lipid[i][0]))
        fl1.write("\t")
        fl1.write(str(min_dist_A_lipid[i][0]))
        fl1.write("\t")
        fl1.write(str(min_dist_B_lipid[i][0]))
        fl1.write("\t")
        fl1.write(str(memb_min_distA_lipid[i][0]))
        fl1.write("\t")
        fl1.write(str(memb_min_distB_lipid[i][0]))
        fl1.write("\t")
        fl1.write(str(rmsd_values_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_BtoA[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(rmsd_values_BtoB[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_B[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(memb_dist[i][0]))
        fl1.write("\t")
        fl1.write(str(Z_angle_A[i][0]))
        fl1.write("\t")
        fl1.write(str(Z_angle_B[i][0]))
        fl1.write("\t")
        fl1.write(str(delta_z_A[i][0]))
        fl1.write("\t")
        fl1.write(str(delta_z_B[i][0]))
        fl1.write("\t")
        fl1.write(str(lipid_count[i][0]))
        fl1.write("\t")
        fl1.write(str(end_to_end_A[i][0]))
        fl1.write("\t")
        fl1.write(str(end_to_end_B[i][0]))
        # fl1.write("\t")
        # fl1.write(str(wall_z_A[i][0]))
        # fl1.write("\t")
        # fl1.write(str(wall_z_B[i][0]))
        fl1.write("\n")

with open("./Delta_Z_measurements.txt",'a') as fl1:
    fl1.write("#Frame\tTimestep\tcomA\tcomB\tdelta_z_A\tdelta_z_B\tz_A_membres\tz_B_membres\tdz_xy_Aupp\tdz_xy_Alow\tdz_xy_Bupp\tdz_xy_Blow\n")
    for i in range(n_frames):
        fl1.write(str(i))
        fl1.write("\t")
        fl1.write(str(u_copy.trajectory[i].time))
        fl1.write("\t")
        fl1.write(str(com_z_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_B[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(delta_z_A[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(delta_z_B[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(com_z_A_hi_flx[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(com_z_B_hi_flx[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(com_z_A_low_flx[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(com_z_B_low_flx[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(com_z_A_med_flx[i][0]/10.0))
        # fl1.write("\t")
        # fl1.write(str(com_z_B_med_flx[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_A_membres[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(com_z_B_membres[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(dz_xy_A_upp[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(dz_xy_A_low[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(dz_xy_B_upp[i][0]/10.0))
        fl1.write("\t")
        fl1.write(str(dz_xy_B_low[i][0]/10.0))
        fl1.write("\n")

with open("./Pitch_Roll_measurements.txt", 'a') as fl2:
    fl2.write("#Frame\tTimestep\tZ_angle_A\tPitchA1\tPitchA2\tPitchA3\tPitchA4\tZ_angle_B\tPitchB1\tPitchB2\tPitchB3\tPitchB4\tRollA1\tRollA2\tRollA3\tRollA4\tRollA5\tRollB1\tRollB2\tRollB3\tRollB4\tRollB5\n")
    for i in range(n_frames):
        fl2.write(str(i))
        fl2.write("\t")
        fl2.write(str(u_copy.trajectory[i].time))
        fl2.write("\t")
        fl2.write(str(Z_angle_A[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A1[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A2[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A3[i][0]))
        fl2.write("\t")
        fl2.write(str(pitch_A4[i][0]))
        fl2.write("\t")
        fl2.write(str(Z_angle_B[i][0]))
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





curr_data.totalFrames = curr_data.totalFrames + n_frames
#Storing the analysis in pickle file
with open(pick_path,"wb") as pick_handle2:
    pickle.dump(curr_data,pick_handle2)

print("Trajectory details: ")
print("Total Frames: ",n_frames)
print("Minimun Distance: ",min(min_dist_lipid))
print("Total frames bound :",curr_data.bound_frames,curr_data.bound_frames/curr_data.totalFrames)
print("Bound Intervals: ",curr_data.boundIntervals)
print(curr_data.inTime,curr_data.outTime)
print(comProA)


