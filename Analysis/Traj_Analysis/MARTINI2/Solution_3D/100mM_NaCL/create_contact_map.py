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
print("No. of frames in Time Dict: ", len(time_vs_con_dict.keys()))
# sys.exit()
cluster_states = list(set(time_dict.values()))
print("Number of Clusters: ",cluster_states)
# print(list(time_vs_con_dict.keys())[-10:])
# sys.exit()
times = list(time_vs_con_dict.keys())
# print(times[679])

time2 = list(time_dict.keys())
# print(time2[0])
# sys.exit()

arr_size = list(time_vs_con_dict.values())[0].shape
clust_contacts = {st : np.zeros(arr_size) for st in cluster_states}
clust_count = {st:0 for st in cluster_states}
time_error = []
"""
"""
max_contact = {}
contact_dict = {}
total_patches_contact = {}
for time,contact_patches in time_vs_con_dict.items():
    try:

        count=0
        st = time_dict[time]
        # print("CLuster: ", st)
        clust_contacts[st]+=contact_patches
    except KeyError:
        time_error.append(time)

print("Unmatched times: ",len(time_error))
# sys.exit()
for st,contact_patches in clust_contacts.items():
    #Normalize contacts
    sum_patch = np.sum(contact_patches,axis=None)
    contact_patches = contact_patches/sum_patch

    fig,ax = plt.subplots()
    hm = ax.imshow(contact_patches)
    ax.set_xticks(np.arange(11))
    ax.set_yticks(np.arange(11))
    ax.set_xlabel("Chain B")
    ax.set_ylabel("Chain A")
    fig.colorbar(hm,ax=ax,aspect=40,label="Fraction of Contacts",orientation='horizontal',panchor=(0.5,0.0),pad=0.2,shrink=0.5)
    ax.set_yticklabels(['d1','d2','d3','M','NP-1','NP-2','NNP-1','NNP-2','tail','NNP-3','NNP-4'])
    ax.set_xticklabels(['d1','d2','d3','M','NP-1','NP-2','NNP-1','NNP-2','tail','NNP-3','NNP-4'])
    plt.setp(ax.get_xticklabels(),rotation=45, ha='center',fontsize=10,va='top')
    plt.setp(ax.get_yticklabels(),va='center',fontsize=10)
    title_str = "Contact Map: " + str(st)
    ax.set_title(title_str)
    plt.show()

# fig,ax = plt.subplots()
# #fig.figsize = [50,50]
# residues_types=key_map_dict.keys()
# hm = ax.imshow(final_contact_array)
# #for i in range(len(key_map_dict.keys())):
# #    residues_list[i] = str(i+1) + residues_list[i]
# ax.set_xticks(np.arange(len(residues_types)))
# ax.set_yticks(np.arange(len(residues_types)))
# fig.colorbar(hm)
# ax.set_xticklabels(residues_types)
# ax.set_yticklabels(residues_types)
# ax.set_xlabel('Chain B')
# ax.set_ylabel('Chain A')
# print(len(residues_list))
# plt.setp(ax.get_xticklabels(),rotation=45, ha='center', rotation_mode='anchor',fontsize=10)
# plt.setp(ax.get_yticklabels(),va='center',fontsize=10)
#
# ax.set_title("Contact Map")
#print(contacts_patches)
# print(max_contact)
# print(contact_dict)
# print(total_patches_contact)

# pick_path2 = "./time_vs_distarray.pickle"
# with open(pick_path2,"rb") as pick_handle2:
#     time_vs_darray = pickle.load(pick_handle2)
#
# time_arr = []
# min_dist = []
# max_dist  =[]
# for time,dist in time_vs_darray.items():
#     time_arr.append(time)
#     min_dist.append(np.amin(dist))
#     max_dist.append(np.amax(dist))
# print(dist.shape)
#
# fig,ax = plt.subplots()
# ax.plot(time_arr,min_dist,linewidth=0.6,label="min")
# ax.plot(time_arr,max_dist,linewidth=0.6,label="max")
# ax.legend()
# ax.set_xlabel("Time")
# ax.set_ylabel("Ang")
# plt.show()
