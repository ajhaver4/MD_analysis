import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import math

from Patches_Class import PatchesClass




"""
Begin Analysis
"""
#The pickle file will store contacts data from each batch of trajectory
import pickle
import os
pick_path = sys.argv[3]   #Path to the pickle file containing contact data
with open(pick_path,"rb") as pick_handle:
    clust_contacts = pickle.load(pick_handle)



label_colours = {0:'forestgreen', 1:'teal', 2:'crimson', 3:'gold', 4:'orchid', 5:'peru', 6:'mediumpurple', 7:'darkorange',8:'steelblue',-1:'olivedrab'}


def convert_to_hist(contact_array):
    contact_count=[]
    for i in range(len(contact_array)):
        # print(contact_array[i])
        if contact_array[i]>=1:
            residue_count = np.ones(int(contact_array[i]))*i
            contact_count+= list(residue_count)
    
    return(contact_count)


patchclass = PatchesClass(sys.argv[1],sys.argv[2])


def patches_map(contact_final,patches_A_dict,patches_B_dict,res_num):
    
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
        for k in range(contact_final[j].shape[0]):
            pB = patches_B_dict[k]
            curr_contact = contact_final[j][k]
            #Increasing count of a contact patch pair corresponding to the residue pair
            contacts_patches[pA-1,pB-1]+= curr_contact

    
    figp,axp = plt.subplots()
    hm = axp.imshow(contacts_patches)
    axp.set_xticks(np.arange(11))
    axp.set_yticks(np.arange(11))
    axp.set_ylabel("Chain A")
    axp.set_xlabel("Chain B")
    fig.colorbar(hm,ax=axp,aspect=40,label="Fraction of Contacts",orientation='horizontal',panchor=(0.5,0.0),pad=0.2,shrink=0.5)
    axp.set_yticklabels(['d1','d2','d3','M','NP-1','NP-2','NNP-1','NNP-2','tail','NNP-3','NNP-4'])
    axp.set_xticklabels(['d1','d2','d3','M','NP-1','NP-2','NNP-1','NNP-2','tail','NNP-3','NNP-4'])
    plt.setp(axp.get_xticklabels(),rotation=45, ha='center',fontsize=10,va='top')
    plt.setp(axp.get_yticklabels(),va='center',fontsize=10)
    title_str = "Contact Map: " + "C-" + str(int(st))
    axp.set_title(title_str)

    



    return(contacts_patches)

def residue_map(contact_final,patch_class):
    
    num_res = patch_class.num_sc
    
    figp,axp = plt.subplots()
    hm = axp.imshow(contact_final)
    axp.set_xticks(np.arange(num_res))
    axp.set_yticks(np.arange(num_res))
    axp.set_ylabel("Chain A")
    axp.set_xlabel("Chain B")
    fig.colorbar(hm,ax=axp,aspect=40,label="Fraction of Contacts",orientation='horizontal',panchor=(0.5,0.0),pad=0.2,shrink=0.5)
    # axp.set_yticklabels(['d1','d2','d3','M','NP-1','NP-2','NNP-1','NNP-2','tail','NNP-3','NNP-4'])
    # axp.set_xticklabels(['d1','d2','d3','M','NP-1','NP-2','NNP-1','NNP-2','tail','NNP-3','NNP-4'])
    plt.setp(axp.get_xticklabels(),rotation=45, ha='center',fontsize=10,va='top')
    plt.setp(axp.get_yticklabels(),va='center',fontsize=10)
    title_str = "Contact Map: " + "C-" + str(int(st))
    axp.set_title(title_str)
    plt.show()


for st,contact_patches in clust_contacts.items():
    fig,[ax1,ax2] = plt.subplots(2,1,figsize=(6,14))
    #Normalize contacts
    sum_patch = np.sum(contact_patches,axis=None)
    contact_patches = contact_patches #/sum_patch
    print("Total Sum of contact: ",sum_patch)

    #Getting 10 most in contact residues
    # flatten_array = contact_patches.flatten()
    # sorted_indx = np.argsort(flatten_array)
    # max_residue_indx = sorted_indx[-10:]
    # max_residue = np.unravel_index(max_residue_indx,contact_patches.shape)


    # patch_contacts = patches_map(contact_patches,patchclass.patches_A_dict,patchclass.patches_B_dict,patchclass.num_sc)
    residue_map(contact_patches,patchclass)
    continue
   

    chainA_sumcontacts = contact_patches.sum(axis=1)
    chainB_sumcontacts = contact_patches.sum(axis=0)
    # for i in range(contact_patches.shape[0]):
    #     chainB_sumcontacts= []
    #     for j in range(contact_patches.shape[1]):
    #         sum_contacts = np.sum(contact_patches[:,j])
    #         chainB_sumcontacts.append(sum_contacts)
        
    residue_histA = convert_to_hist(chainA_sumcontacts)
    residue_histB = convert_to_hist(chainB_sumcontacts)

    # print("Residue Hist A: ",residue_histA)

        
    # print("Count from Residues Hist: ", residue_histA.count(104))
    # print("COutn from Direct Sum: ", chainA_sumcontacts[104])

    

    ###Background Coloring according to Patch Classification
    # for patch,resids in patchclass.patch_labels.items():
    #     print("Patch: ",patch)
    #     print("Residues: ",resids)

    #     if patch == 'd1' or patch =='d2' or patch == 'd3' or patch == 'NP-1' or patch == 'NP-2':
    #         bg_color='magenta'
    #     elif patch == 'M':
    #         bg_color='black'
    #     elif patch == 'NNP-1' or patch == 'NNP-2' or patch == 'NNP-3' or patch == 'NNP-4':
    #         bg_color='cyan'
    #     elif patch == 'tail':
    #         bg_color='yellow'
    #     for res in resids:
            
    #         ax.axvspan(res-1,res,alpha=0.2,color=bg_color)

    #### Background coloring according to Biochemical Classification
    count=0

    biochem_sum_count = {key:[0,0] for key in patchclass.biochem_dict.keys()}
    # biochem_sum_countB = {key:0 for key in patchclass.biochem_dict.keys()}
    biochem_histA=[]
    biochem_histB=[]
    for patch,resids in patchclass.biochem_dict.items():
        # print("Patch: ",patch)
        # print("Residues: ",resids)

        if patch == 'Polar':
            bg_color='cyan'
        elif patch == '+ Charge':
            bg_color='lime'
        elif patch == '- Charge':
            bg_color='magenta'
        elif patch == 'Non-Polar':
            bg_color='yellow'
        
        # biochem_histA=residue_histA
        # biochem_histB=residue_histB

        

        # for res in resids:
            
        #     ax.axvspan(count,count+1,alpha=0.2,color=bg_color)
        #     for r in range(len(residue_histA)):
        #         if res == residue_histA[r]:
        #             biochem_histA[r] = count
        #             biochem_sum_count[patch][0]+=1
            
        #     for r in range(len(residue_histB)):
        #         if res == residue_histB[r]:
        #             biochem_histB[r] = count
        #             biochem_sum_count[patch][1]+=1
            
        #     count+=1

        patch_countA=0
        patch_countB=0
        for res in resids:
            
            ax1.axvspan(count,count+1,alpha=0.2,color=bg_color)
            res_countA = chainA_sumcontacts[res]
            biochem_histA+= [count]*int(res_countA)
            patch_countA+=res_countA

            res_countB = chainB_sumcontacts[res]
            biochem_histB+= [count]*int(res_countB)
            patch_countB+=res_countB
            
            count+=1
        
        biochem_sum_count[patch][0] = patch_countA
        biochem_sum_count[patch][1] = patch_countB
        print("Total residues counted: ",count)
        print("Size of Contact Array : ",np.shape(contact_patches))
    
    # print("Modified Residue Hist A: ",biochem_histA)
    # residue_histA = biochem_histA
    # residue_histB = biochem_histB

    #Plotting in sorted order of biochem classification
    bin_width = np.arange(0,256,1)
    plot_bioA=ax1.hist(biochem_histA,bins=bin_width,color=label_colours[int(st)],alpha=0.8,label='C-' + str(int(st)))
    ax12=ax1.twinx()
    plot_bioB=ax12.hist(biochem_histB,bins=bin_width,color=label_colours[int(st)],alpha=0.8)
    ax12.invert_yaxis()
    ax1.legend()

    ax1.set_xlabel("Residue Number")
    ax1.set_ylabel("Fraction of Contacts")
    title_lb = "Contact Histogram: " + "C-" + str(int(st))
    ax1.set_title(title_lb)

    #Plotting Normal residue level contact
    bin_width = np.arange(0,256,1)
    plotA=ax2.hist(residue_histA,bins=bin_width,color=label_colours[int(st)],alpha=0.8,label='C-' + str(int(st)))
    ax22=ax2.twinx()
    plotB=ax22.hist(residue_histB,bins=bin_width,color=label_colours[int(st)],alpha=0.8)
    ax22.invert_yaxis()
    ax2.legend()

    ax1.set_xlabel("Residue Number")
    ax1.set_ylabel("Fraction of Contacts")
    title_lb = "Contact Histogram: " + "C-" + str(int(st))
    ax1.set_title(title_lb)


    color_dict = {'Polar':'cyan','+ Charge':'lime','- Charge':'magenta','Non-Polar':'yellow'}
    fig2,ax2=plt.subplots()
    for patch,count in biochem_sum_count.items():
        print(count)
        ax2.bar(patch,count[0],color=color_dict[patch],label=patch,alpha=0.8)
        ax2.bar(patch,count[1],color=color_dict[patch],label=patch,bottom=count[0],alpha=0.4,hatch='/')
    
    # for patch,count in biochem_sum_countB.items():
    #     ax2.bar(patch,count,color=color_dict[patch],label=patch)

    plt.show()


