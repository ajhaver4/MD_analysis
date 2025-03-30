# import pandas
import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
# from sklearn.cluster import KMeans
from matplotlib.markers import MarkerStyle as markerstyle
from mpl_toolkits.mplot3d import axes3d
from MDAnalysis.lib.distances import distance_array as dist_array
from MDAnalysis.analysis import contacts
import time as time_mod
"""
Input Information IMPORTANT!!!!!!!
#First input has to a .tpr file. Since numbering and residue types are different if .gro file is used!!!!!!!!
In .tpr file numbering starts from 0 and .gro file numbering starts from 1.
"""
#clust_file = sys.argv[1]
struct_file = sys.argv[1]
traj_file = sys.argv[2]

#Define MDAnalysis Universe
#u = mda.Universe(struct_file,traj_file)
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()

protein = u.select_atoms("byres name BB")
protein_rms = u.select_atoms("byres name BB")
# print(protein)

n = len(protein)
n_frames = u.trajectory.n_frames

chainA = protein.select_atoms("bynum 1-469")
chainB = protein.select_atoms("bynum 470-938")
chainB_mic = protein_rms.select_atoms("bynum 470-938")
print(chainA)
sys.exit()

A_sc_selection = "(not name BB) and bynum 1-469"
B_sc_selection = "(not name BB) and bynum 470-938"

chainA_sc = chainA.select_atoms(A_sc_selection)
chainB_sc = chainB_mic.select_atoms(B_sc_selection)
# print(chainA_sc)

num_sc = len(chainA_sc)
num_total = len(chainA)
print(num_sc)



"""
CONTACT PATCH DEFINITIONS
"""
#Defining Patches:
d1_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 0:4 or resid 6:11 or resid 13)")
memb_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 5 or resid 12 or resid 15 or resid 19 or resid 75 or resid 79 or resid 82)")
patch1_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 14 or resid 16:18 or resid 20:29)")
patch2_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 30:54)")
patch3_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 55:74 or resid 76:78 or resid 80:81 or resid 83:99)")
tail_A=protein.select_atoms("(not name BB) and bynum 1-469 and (resid 100:135)")
NN9_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 136:154)")
d3_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 155:162)")
patch11_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 163:173)")
d2_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 174:178)")
patch13_A = protein.select_atoms("(not name BB) and bynum 1-469 and (resid 179:213)")

patch1_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 214:218 or resid 220:225 or resid 227:228 or resid 230:232 or resid 234:240)")
memb_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 219 or resid 226 or resid 229 or resid 233 or resid 289 or resid 293 or resid 296)")
d1_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 241:255)")
patch2_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 256:268)")
patch3_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 269:288 or resid 290:292 or resid 294:295 or resid 297:313)")
tail_B=protein.select_atoms("(not name BB) and bynum 470-938 and (resid 314:351)")
d2_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 352:359)")
patch9_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 360:368)")
d3_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 369:376)")
NN11_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 377:405)")
patch13_B = protein.select_atoms("(not name BB) and bynum 470-938 and (resid 406:427)")

patches_A_dict= {}
patches_B_dict = {}
count=0
for i in range(len(chainA_sc)):
    res = chainA_sc[i].resid

    if res in d1_A.resids:
        patches_A_dict[i]=1
    elif res in d2_A.resids:
        patches_A_dict[i]=2
    elif res in d3_A.resids:
        patches_A_dict[i]=3
    elif res in memb_A.resids:
        patches_A_dict[i]=4
    elif res in NN9_A.resids:
        patches_A_dict[i]=5
    elif res in patch2_A.resids:
        patches_A_dict[i]=6
    elif res in patch1_A.resids:
        patches_A_dict[i]=7
    elif res in patch3_A.resids:
        patches_A_dict[i]=8
    elif res in tail_A.resids:
        patches_A_dict[i]=9
    elif res in patch11_A.resids:
        patches_A_dict[i]=10
    elif res in patch13_A.resids:
        patches_A_dict[i]=11

for i in range(len(chainB_sc)):
    res = chainB_sc[i].resid
    if res in d1_B.resids:
        patches_B_dict[i]=1
    elif res in d2_B.resids:
        patches_B_dict[i]=2
    elif res in d3_B.resids:
        patches_B_dict[i]=3
    elif res in memb_B.resids:
        patches_B_dict[i]=4
    elif res in NN11_B.resids:
        patches_B_dict[i]=5
    elif res in patch1_B.resids:
        patches_B_dict[i]=6
    elif res in patch2_B.resids:
        patches_B_dict[i]=7
    elif res in patch3_B.resids:
        patches_B_dict[i]=8
    elif res in tail_B.resids:
        patches_B_dict[i]=9
    elif res in patch9_B.resids:
        patches_B_dict[i]=10
    elif res in patch13_B.resids:
        patches_B_dict[i]=11

print(patches_A_dict)
print(patches_B_dict)
print("LENGTH OF PATCHESA: ", len(patches_A_dict.keys()))
print("LENGTH OF PATCHESB: ", len(patches_B_dict.keys()))
sys.exit()
#num_res=4
"""
-------------------------------------------------------------------------------
---------------------------- METHOD DEFINITIONS -------------------------------
-------------------------------------------------------------------------------
"""

class ContactAtom():
    def __init__(self,atom,res_len):
        self.atom = atom
        self.contact_array = np.zeros((1,res_len))

    def add_contact(self,con_array):
        self.contact_array = self.contact_array + con_array

    #def __init__(self,contact_array):
        #self.contact_array = contact_array
def create_contact_map(final_dist_array,con_radius,res_num,patches_A_dict,patches_B_dict):
    global global_frames
    global curr_clust_contacts
    contact_final = np.zeros((res_num,res_num))
    contacts_patches = np.zeros((len(set(patches_A_dict.values())),len(set(patches_B_dict.values()))))
    #con_radius = 10
    #for i in range(len(curr_clust_contacts.keys())):
        #print("Cluster: ",i)
        #curr_clust_contacts[i][0,:] = np.ones((1,num_res))*(1)
        #print(curr_clust_contacts[i])
    #Create residue wise contact map for whole trajectory
    # print("Cluster STATE: ",state)
    #print(np.max(final_dist_array))
    res = 0
    clust_list=[]
    for j in range(res_num):
        #Loop over 1st dim i.e. chainA residues
        pA = patches_A_dict[j]
        res=res+1
        #Calculate contact matrix
        con = contacts.contact_matrix(final_dist_array[j],radius=con_radius)
        contact_final[j][:]=con.astype(float)
        print(len(contact_final[j].shape[0]))
        for k in range(contact_final[j].shape[0]):
            pB = patches_B_dict[k]
            curr_contact = contact_final[j][k]
            #Increasing count of a contact patch pair corresponding to the residue pair
            contacts_patches[pA-1,pB-1]+= curr_contact

    return(contacts_patches)

def check_image(pos1,pos2):

    """
    This function checks if the current image of chain B is the minimum distance image from chain A or not.
    It checks by comparing the distance between com_A and com_B in each dimension.
    If it is greater than half_box len then it creates a translation vector that shifts chain B along that dimension in the opposite direction.
    As explained in the code below, this is required because while calc. RMSD the distances have to be calculated using the Minimum Image of chain B

    Returns a lot values but the important ones are the translation vec and the boolean mic which says that chain B is being translated.
    """
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
    """
    Calculates the relative position vector w.r.t to COM
    """
    pos_array = Ag1.positions
    return(pos_array - com)

final2 = np.zeros((num_sc,num_sc))
ts = u.trajectory.ts
global box_dims
box_dims = ts.dimensions

"""
-------------------------------------------------------------------------------
---------------------------- CREATING PICKLE FILE/LOADING DATA FROM PICKLE -----
-------------------------------------------------------------------------------
"""
#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
pick_path = "./time_vs_patches.pickle"

try:
    if not os.path.exists(pick_path):
        time_vs_con_dict = {}
    else:
        with open(pick_path,"rb") as pick_handle:
            time_vs_con_dict = pickle.load(pick_handle)
except:
    pick_path = "./time_vs_patches_" + str(ts.time) + "_us.pickle"

# pick_path2 = "./time_vs_distarray.pickle"
# if not os.path.exists(pick_path2):
#     time_vs_darray = {}
# else:
#     with open(pick_path2,"rb") as pick_handle2:
#         time_vs_darray = pickle.load(pick_handle2)

time_vs_darray = {}
"""
-------------------------------------------------------------------------------
---------------------------- LOOP OVER TRAJECTORY -----------------------------
-------------------------------------------------------------------------------
"""


print("Reading trajectory frames...")
print("Total frames: ", n_frames)
time_error=[]
#print(clust_contacts.keys())
#num_res=4
#print(clust_contacts)


step = 20
chunks = np.arange(0,n_frames,step)
fr_count = 0
count=0
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
        # ts = u.trajectory.ts
        protein.unwrap(reference='cog',inplace=True)
        # protein_rms.unwrap(reference='cog',inplace=True)
        #print(ts)
        ts = u.trajectory.ts
        try:
            comProA = chainA.centroid(pbc=False,compound='segments')
            comProB = chainB.centroid(pbc=False,compound='segments')
            (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
            relVec = calc_relVec(chainB,comProB)
            new_comProB = comProB + transVec
            new_positions_B = relVec + new_comProB
            new_array = ts.positions
            new_array[469:938,:] = new_positions_B

            if mic:
                print("********   MIC  **********")
                # print("Frame: ",fr_count)
                # print("Before: ", protein.positions[547,:])
                # print("Rel Vector: ",new_positions_B[0,:])
                #print("Timestep: ", ts.positions[458,:])
                # protein_rms.positions[547:1094,:] = new_positions_B
                chainB_mic.translate(transVec)
                # print("After Protein_rms: ", protein.positions[547,:])
                # print("After chainB_mic: ", chainB.positions[0,:])

            print("Timestep: ",ts.time)
            print("Copy Timestep: ",u_copy.trajectory[fr_count].time)
            darray = dist_array(chainA.positions,chainB_mic.positions,box=box_dims)

            print(np.amin(darray))

            if np.amin(darray) < 8.0:
                print("Contact")
                sc_dist_array = dist_array(chainA_sc.positions,chainB_sc.positions,box=box_dims)
                final2[:][:] = sc_dist_array
                patches_contacts_array = create_contact_map(final2,6,num_sc,patches_A_dict,patches_B_dict)
                time_vs_con_dict[u_copy.trajectory[fr_count].time]=patches_contacts_array
                time_vs_darray[u_copy.trajectory[fr_count].time]=darray
                count+=1
        except KeyError:
            time_error.append(ts)
            # fr_count+=1
            print("ERROR: ",time_error)
        except Exception as e:
            # fr_count+=1
            print("There was a error: ", e)

        # if i>10:
        #     for i in range(len(curr_clust_contacts.keys())):
        #         print("Cluster: ",i)
        #         #curr_clust_contacts[i][0,:] = np.ones((1,num_res))*(1)
        #         print(np.max(curr_clust_contacts[i]))
        #     sys.exit()
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
print("Number of frames in contact: %.2f " %(count*100/n_frames))

print(time_vs_con_dict.keys())

with open(pick_path,"wb") as write_handle:
    pickle.dump(time_vs_con_dict,write_handle)

time_str = int(u_copy.trajectory[fr_count].time/1e6)
pick_path2 = "./time_vs_distarray" + str(time_str) + "_us.pickle"
with open(pick_path2,"wb") as write_handle2:
    pickle.dump(time_vs_darray,write_handle2)
