# import pandas
import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
# from sklearn.cluster import KMeans
from matplotlib.markers import MarkerStyle as markerstyle
from mpl_toolkits.mplot3d import axes3d
from Clust_times import time_dict
from MDAnalysis.lib.distances import distance_array as dist_array
from MDAnalysis.analysis import contacts
import time as time_mod

#clust_file = sys.argv[1]
struct_file = sys.argv[1]
traj_file = sys.argv[2]

#Define MDAnalysis Universe
#u = mda.Universe(struct_file,traj_file)
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()

protein = u.select_atoms("byres name BB")

n = len(protein)
n_frames = u.trajectory.n_frames

chainA = protein.select_atoms("bynum 1-469")
chainB = protein.select_atoms("bynum 470-938")

#print(chain_A.select_atoms("resname ARG ")[1:4])
#print(chain_B)

chainA_bb = chainA.select_atoms("type B")
chainB_bb = chainB.select_atoms("type B")

A_sc_selection = "(byres name BB) and not name BB"
B_sc_selection = "(byres name BB) and not name BB"

chainA_sc = chainA.select_atoms(A_sc_selection)
chainB_sc = chainB.select_atoms(B_sc_selection)

global num_res
num_res = len(chainA_sc)
print(num_res)

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
def create_contact_map(final_dist_array,con_radius,res_num,state):
    global global_frames
    global num_res
    global curr_clust_contacts
    contact_final = np.zeros((res_num,res_num))
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
        res=res+1
        con = contacts.contact_matrix(final_dist_array[j],radius=con_radius)
        contact_final[j][:]=con.astype(float)
        #print(con.astype(float))
        curr_clust_contacts[state][j,:] = curr_clust_contacts[state][j,:] + contact_final[j][:].reshape((1,res_num))
    #curr_clust_contacts[0][0,:] = np.ones((1,num_res))*(1)
        #clust_list.append(state)
    # print("No. of residues: ",res)
    #print("Clusters: ",set(clust_list))
    #print(curr_clust_contacts

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

"""
-------------------------------------------------------------------------------
---------------------------- CREATING PICKLE FILE/LOADING DATA FROM PICKLE -----
-------------------------------------------------------------------------------
"""
#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
pick_path = "./clust_contacts_data.pickle"
cluster_states = list(set(time_dict.values()))
if not os.path.exists(pick_path):
    contact_array = np.zeros((num_res,num_res))
    clust_contacts = {st : np.zeros((num_res,num_res)) for st in cluster_states}
else:
    with open(pick_path,"rb") as pick_handle:
        clust_contacts = pickle.load(pick_handle)

cluster_states = list(set(time_dict.values()))
contact_array = np.zeros((num_res,num_res))
#clust_contacts = {st : np.zeros((num_res,num_res)) for st in cluster_states}
global curr_clust_contacts
curr_clust_contacts = clust_contacts
clust_count = {st:0 for st in cluster_states}

"""
-------------------------------------------------------------------------------
---------------------------- LOOP OVER TRAJECTORY -----------------------------
-------------------------------------------------------------------------------
"""

final2 = np.zeros((num_res,num_res))
ts = u.trajectory.ts
global box_dims
box_dims = ts.dimensions
print("Reading trajectory frames...")
time_error=[]
#print(clust_contacts.keys())
#num_res=4
#print(clust_contacts)


step = 200
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
                # chainB.translate(transVec)
                # print("After Protein_rms: ", protein.positions[547,:])
                # print("After chainB_mic: ", chainB.positions[0,:])

            print("Timestep: ",ts.time)
            cluster = time_dict[ts.time]
            clust_count[cluster]+=1
            d_array_2 = dist_array(chainA_sc.positions,chainB_sc.positions,box=box_dims)
            final2[:][:] = d_array_2
            #print(np.max(d_array_2))
            create_contact_map(final2,6,num_res,cluster)
            count+=1
            #print("Frame: ",i)
            #print("Contact Array")
            #print(np.max(curr_clust_contacts[0]))
            #print(clust_contacts)
            # fr_count+=1
        except KeyError:
            time_error.append(ts)
            # fr_count+=1
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


print("******------------------*******")
print("******------------------*******")
for i in range(len(curr_clust_contacts.keys())):
    print("Cluster: ",i)
    #curr_clust_contacts[i][0,:] = np.ones((1,num_res))*(1)
    print(curr_clust_contacts[i])
    print("Max Contact: ", np.max(curr_clust_contacts[i]))
    print("No. of frames: ",clust_count[i])
print("No. of frames which gave a time error: ", len(time_error))
print("Count: ",count)
#sys.exit()
print("******------------------*******")
print("******------------------*******")


#Contact Map
print("Creating Contact map")
max_res_dict = {st : [] for st in cluster_states}
print("Cluster Contact Data")
print(np.max(curr_clust_contacts[0]))
print("CLuster 2")
print(np.max(curr_clust_contacts[2]))
for clust,contact_list in curr_clust_contacts.items():
    max_contact_res = ""
    max_contact = 0
    for i in range(len(contact_list)):
        res_net_contacts = np.sum(contact_list[i,:],axis=0)
        #print("Shape: ", contact_list[i,:].shape)
        curr_contact_res = np.argmax(contact_list[i,:])
        #print(curr_contact_res)
        curr_contact = contact_list[i,curr_contact_res]
        if curr_contact > max_contact:
            max_contact=curr_contact

            res_A = str(chainA_sc[i].resname) + str(i)
            res_B = str(chainB_sc[curr_contact_res].resname) + str(curr_contact_res)
            max_contact_res=(res_A,res_B)
        # for j in range(len(res_net_contacts)):
        #     final_contact_array[key_map_dict[res_type_dict[i]]][key_map_dict[res_type_dict[j]]] =\
        #      final_contact_array[key_map_dict[res_type_dict[i]]][key_map_dict[res_type_dict[j]]]+res_net_contacts[j]
    print(clust,max_contact_res,max_contact)
    max_res_dict[clust] = (max_contact_res,max_contact)
    print("Max:",np.max(contact_list))
#Storing the analysis in pickle file
with open(pick_path,"wb") as pick_handle2:
    pickle.dump(curr_clust_contacts,pick_handle2)

print(curr_clust_contacts[0].shape)
print(max_res_dict)
