import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import math
from Clust_times import time_dict



"""
Begin Analysis
"""
#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
pick_path = "./time_vs_patches.pickle"
with open(pick_path,"rb") as pick_handle:
    time_vs_con_dict = pickle.load(pick_handle)

print("No. of frames in Contact Dict: ", len(time_vs_con_dict.keys()))
print("No. of frames in Time Dict: ", len(time_dict.keys()))
# sys.exit()
cluster_states = list(set(time_dict.values()))
print("Number of Clusters: ",cluster_states)
# print(list(time_vs_con_dict.keys())[-10:])
# sys.exit()
times = list(time_vs_con_dict.keys())
# print(times[679])

time2 = list(time_dict.keys())


def get_value_from_dict(time,time_vs_con_dict):
    all_times = np.array(list(time_vs_con_dict.keys()))
    print(all_times)
    diff = all_times-time
    diff = np.abs(diff)
    if np.min(diff) > 200:
        print("Time not found in dictionary")
        return(None)
    min_indx = np.argmin(diff)
    actual_time = all_times[min_indx]
    print(actual_time,time)
    return(time_vs_con_dict[actual_time])

arr_size = list(time_vs_con_dict.values())[0].shape

clust_count = {st:0 for st in cluster_states}
time_error = []

#Checking if previous contacts are stored in a pickle file
print("Checking if previous contacts are stored in a pickle file...")
print(len(sys.argv))
if len(sys.argv) > 1:
    with open(sys.argv[1],"rb") as pick_handle:
        clust_contacts = pickle.load(pick_handle)
else:
    print("No pickle file found. Initializing new dictionary...")
    clust_contacts = {st : np.zeros(arr_size) for st in cluster_states}

"""
"""
max_contact = {}
contact_dict = {}
total_patches_contact = {}
contact_frames=0
#for time,contact_patches in time_vs_con_dict.items():
#    try:
#
#        count=0
#
#        st = time_dict[time]
#        # print("CLuster: ", st)
#        clust_contacts[st]+=contact_patches
#    except KeyError:
#        print(time)
#        time_error.append(time)

for time,contact_patches in time_vs_con_dict.items():
    try:
        st = time_dict[time]
        # print("CLuster: ", st)
        clust_contacts[st]+=contact_patches
        contact_frames+=1
         
        try:
            contact_prev = get_value_from_dict(time-200,time_vs_con_dict)
            contact_next = get_value_from_dict(time+200,time_vs_con_dict)
            if contact_prev is None or contact_next is None:
                continue
            total_contacts = contact_prev + contact_next + contact_patches
            total_contacts = total_contacts/3
            clust_contacts[st]+=total_contacts.astype('int')
            contact_frames+=1
        except Exception as e:
            print(e)
            print("Error in getting previous and next contacts")
#            sys.exit()
            continue
    except KeyError:
        print(time)
        time_error.append(time)


        

print("Unmatched times: ",len(time_error))
print("Frames Counted for Contact: ",contact_frames)
# print(time_error)
# sys.exit()


pick_path2 = "./Final_ClusterContacts_Averaged" + ".pickle"
with open(pick_path2,"wb") as write_handle2:
    pickle.dump(clust_contacts,write_handle2)

# label_colours = {0:'forestgreen', 1:'teal', 2:'crimson', 3:'gold', 4:'orchid', 5:'peru', 6:'mediumpurple', 7:'darkorange',8:'steelblue',-1:'olivedrab'}

# for st,contact_patches in clust_contacts.items():
#     #Normalize contacts
#     sum_patch = np.sum(contact_patches,axis=None)
#     contact_patches = contact_patches/sum_patch

#     flatten_array = contact_patches.flatten()
#     sorted_indx = np.argsort(flatten_array)
#     max_residue_indx = sorted_indx[-10:]
#     max_residue = np.unravel_index(max_residue_indx,contact_patches.shape)



#     fig,ax = plt.subplots()
#     for i in range(contact_patches.shape[0]):
#         for j in range(contact_patches.shape[1]):
#             ax.bar(j,np.sum(contact_patches[:,j]))
    
#     plt.show()


