import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import math
from Patches_Class import PatchesClass
from combine_pickles import combine_dict

"""
Begin Analysis
"""

patchclass = PatchesClass(sys.argv[1],sys.argv[2])

#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
if len(sys.argv[4:]) > 1:
    print("Combining multiple pickle files...")
    time_vs_con_dict = combine_dict(sys.argv[3:])
else:

    pick_path = sys.argv[4]
    with open(pick_path,"rb") as pick_handle:
        time_vs_con_dict = pickle.load(pick_handle)



print("No. of frames in Contact Dict: ", len(time_vs_con_dict.keys()))


times = list(time_vs_con_dict.keys())

start_time = np.min(times)
end_time = np.max(times)
print("Time range: %.2f - %.2f us" %(start_time/1e6, end_time/1e6))


def get_value_from_dict(time,time_vs_con_dict):
    all_times = np.array(list(time_vs_con_dict.keys()))
    
    diff = all_times-time
    
    diff = np.abs(diff)
    
    if np.min(diff) > 200:
        print("Time not found in dictionary")
        return
    min_indx = np.argmin(diff)
    
    actual_time = all_times[min_indx]
    
    return(time_vs_con_dict[actual_time])

arr_size = list(time_vs_con_dict.values())[0].shape

#Checking if previous contacts are stored in a pickle file
# print("Checking if previous contacts are stored in a pickle file...")
# if sys.argv[1]:
#     with open(sys.argv[1],"rb") as pick_handle:
#         clust_contacts = pickle.load(pick_handle)
# else:
#     print("No pickle file found. Initializing new dictionary...")
#     clust_contacts = {st : np.zeros(arr_size) for st in cluster_states}

"""
Declaring necessary dictionaries to store contacts informations. 
The initialization depends upon if the analysis is being resumed from a previous chunk of trajectory. 
If the analysis is being resumed then the previous contacts are loaded from the pickle file.
"""



pick_path2 = "./Charged_Interactions_Array" + ".pickle"
pick_path3 = "./Charged_Interactions_Sum" + ".pickle"
pick_path4 = "./Charged_Interactions_Time" + ".pickle"

bool_restart = sys.argv[3]

if bool_restart=='True':
    print("Reading previous contacts from pickle files...")
    if os.path.exists(pick_path2):
        with open(pick_path2,"rb") as pick_handle2:
            contact_dict = pickle.load(pick_handle2)
#        os.remove(pick_path2)
    if os.path.exists(pick_path3):
        with open(pick_path3,"rb") as pick_handle3:
            win_vs_contact = pickle.load(pick_handle3)
            win_keys = list(win_vs_contact.keys())
            #sml_win=win_keys[0]
            med_win=win_keys[0]
            large_win=win_keys[1]
#        os.remove(pick_path3)
    if os.path.exists(pick_path4):
        with open(pick_path4,"rb") as pick_handle4:
            chargepair_timedict = pickle.load(pick_handle4)
#        os.remove(pick_path4)

else:
    contact_dict = {}
    chargepair_timedict={}
    large_win = 1          #In each direction. More like half window size. Total frames covered = 2*window_size+1
    sml_win=0
    med_win=1
    win_vs_contact = {sml_win:{},med_win:{},large_win:{}}

"""
Loop through each timeframe and contact dictionary to look at individual pair-wise contacts
"""
contact_frames=0

def store_contacts(contact_array,window_size):
    #Analyzing each charged-charged interaction
    pos_resids = patchclass.biochem_dict['+ Charge']
    neg_resids = patchclass.biochem_dict['- Charge']
    # print(pos_resids,neg_resids)
    #Chain A residues are accessed by the first dimension of the contact array
    for pos_res in pos_resids:
        for neg_res in neg_resids:
            #Charge pairs where Positive charged residue belongs to A and neg charged residue is on B
            pair_contact = contact_array[pos_res,neg_res]

            #Storing contact information for medium window contacts
            #each timepoint is stored. Kind of expensive but after this I can just access the dictionary(as a pickle) to get the contact information without having to run this.  
            #So the contact vs time information is only stored for one window size i.e. the medium lenght one

            #If there is contact between the pair of residues then add it to the final summed contacts
            if pair_contact >= 1:
                if (pos_res,neg_res) not in win_vs_contact[window_size].keys():
                    win_vs_contact[window_size][(pos_res,neg_res)] = 0
                win_vs_contact[window_size][(pos_res,neg_res)]+=1

            #Charge pairs where Positive charged residue belongs to B and neg charged residue is on A
            contact_other = contact_array[neg_res,pos_res]
                       
            if contact_other >= 1:
                if (neg_res,pos_res) not in win_vs_contact[window_size].keys():
                    win_vs_contact[window_size][(neg_res,pos_res)] = 0
                win_vs_contact[window_size][(neg_res,pos_res)]+=1

print("Largest Window Size: ",large_win)
def store_timemap(curr_time,contact_array,windown_size):
    #Analyzing each charged-charged interaction
    pos_resids = patchclass.biochem_dict['+ Charge']
    neg_resids = patchclass.biochem_dict['- Charge']
    # print(pos_resids,neg_resids)
    #Chain A residues are accessed by the first dimension of the contact array
    for pos_res in pos_resids:
        for neg_res in neg_resids:
            #Charge pairs where Positive charged residue belongs to A and neg charged residue is on B
            pair_contact = contact_array[pos_res,neg_res]
            #Storing contact information for medium window contacts
            #each timepoint is stored. Kind of expensive but after this I can just access the dictionary(as a pickle) to get the contact information without having to run this.  
            #So the contact vs time information is only stored for one window size i.e. the medium lenght one
            if pair_contact >= 1:       #Only store information ifthe residue pair in contact for atleast this medium window size 
                if (pos_res,neg_res) not in contact_dict.keys():
                    contact_dict[(pos_res,neg_res)] = []
                contact_dict[(pos_res,neg_res)].append(pair_contact)

                if curr_time not in chargepair_timedict.keys():
                    chargepair_timedict[curr_time] = []
                chargepair_timedict[curr_time].append((pos_res,neg_res))

            #Charge pairs where Positive charged residue belongs to B and neg charged residue is on A
            contact_other = contact_array[neg_res,pos_res]
                       
            if contact_other >= 1:
                if (neg_res,pos_res) not in contact_dict.keys():
                    contact_dict[(neg_res,pos_res)] = []
                contact_dict[(neg_res,pos_res)].append(contact_other)
                
                if curr_time not in chargepair_timedict.keys():
                    chargepair_timedict[curr_time] = []
                chargepair_timedict[curr_time].append((neg_res,pos_res))


for time,contact_patches in time_vs_con_dict.items():
    
    try:
        #Instead of each frame controbuting to the contact, we will only consider contacts which remain for atleast 3 frames or 600ps. 
        # To do this, take average of 3 frames.(Running average with a window size of 3 frames) and only consider contact if avg >=1.
        #Add adjustable window size for this.
        print("Initializing all contacts")
        contacts_curr = 0    
        contacts_curr = 0
        total_contacts = contacts_curr
        print("Current Time: ",time,np.max(contacts_curr))
        contact_win_sml = np.zeros(arr_size)
        contact_win_med = np.zeros(arr_size)
        contact_prev = 0
        contact_next = 0
        if np.sum(contact_patches) ==0:
            print("No contacts found for this frame. Skipping to the next frame : ",time)
            continue
        for i in range(1,large_win+1):
            # print("Window Size: ",i)
            prev_time = time-(i*200)
            next_time = time+(i*200)
            contact_prev = get_value_from_dict(prev_time,time_vs_con_dict)
            contact_next = get_value_from_dict(next_time,time_vs_con_dict)

            #print("Current Contact: ",contacts_curr[:10,:10])
            #print("Prev Contact Time : ",prev_time,contact_prev[:10,:10])
            #print("Next Contact Time: ",next_time,contact_next[:10,:10])
            
            if contact_prev is None or contact_next is None:
                print("No contacts found for this window size - ",i, "at time: ",time)
                break

            else:

            
                # print("All contacts found for the curent window size. Adding to the total contacts")
                #If prev and next contacts exists then add to the total count. Or else exit the inner loop and start the next iteration of the outer for loop (for-else statement)
                print("Max value in curr frame: ",np.max(contact_patches))
                print("Max value in prev frame: ",np.max(contact_prev))
                print("Max value in next frame: ",np.max(contact_next))
                contacts_curr=contact_patches+contact_prev+contact_next
                #If window reaches a size of 400ps i.e. i=1, then store the contact array
                if i == sml_win:
                    contact_win_sml = contacts_curr/(2*sml_win+1)
                    if np.isinf(contact_win_sml).any():
                        print("Infinite contact found for this window size - ",i, "at time: ",time)
                        print(contact_patches)
                        print(2*sml_win+1)
                        break
                    store_contacts(contact_win_sml,i)
                #If window reaches a size of 2ns i.e. i=5, then store the contact array
                if i == med_win:
                    contact_win_med = contacts_curr/(2*med_win+1)
                    if np.isinf(contact_win_med).any():
                        print("Infinite contact found for this window size - ",i, "at time: ",time)
                        print(contact_patches)
                        print(2*med_win+1)
                        break
                    store_contacts(contact_win_med,i)
                    store_timemap(time,contact_win_med,i) 
            
        
        else:
            #This block is exceuted only if the for loop completes without a break statement
            #It means the largest window size used gathered all contacts for each timepoint without error and can be stored now.
            print("Contacts found for the largest window size - ", i, "at time: ",time)
            contacts_curr = contacts_curr/(2*large_win+1)
            store_contacts(contacts_curr,i)          
            
        
        contact_frames+=1

        if contact_frames%5000==0:
            print("Frames Counted for Contact: ",contact_frames)
            print("Current Time: ",time)
    except Exception as e:
        print("Error in getting previous and next contacts")
        print(e)
        continue



        

# print("Unmatched times: ",len(time_error))
print("Frames Counted for Contact: ",contact_frames)

#Store data in pcikle files
with open(pick_path2,"wb") as write_handle2:
    pickle.dump(contact_dict,write_handle2)


with open(pick_path3,"wb") as write_handle2:
    pickle.dump(win_vs_contact,write_handle2)

with open(pick_path4,"wb") as write_handle2:
    pickle.dump(chargepair_timedict,write_handle2)
