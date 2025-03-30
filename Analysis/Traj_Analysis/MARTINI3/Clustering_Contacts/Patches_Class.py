import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import math
import numpy as np

class PatchesClass:

    def __init__(self,tpr_file,traj_file) -> None:
        """
        Input Information IMPORTANT!!!!!!!
        #First input has to a .tpr file. Since numbering and residue types are different if .gro file is used!!!!!!!!
        In .tpr file numbering starts from 0 and .gro file numbering starts from 1.
        """
        
        
        #Define MDAnalysis Universe
        #u = mda.Universe(struct_file,traj_file)
        self.u = mda.Universe(tpr_file,traj_file)
        self.u_copy = self.u.copy()

        self.protein = self.u.select_atoms("byres name BB")
        self.protein_rms = self.u.select_atoms("byres name BB")

        self.n = len(self.protein)
        
        self.chainA = self.protein.select_atoms("bynum 1-469")
        self.chainB = self.protein.select_atoms("bynum 470-938")
        

        A_sc_selection = "(not name BB) and bynum 1-469"
        B_sc_selection = "(not name BB) and bynum 470-938"

        self.chainA_sc = self.chainA.select_atoms(A_sc_selection)
        self.chainB_sc = self.chainB.select_atoms(B_sc_selection)
        print(self.chainA_sc)
        print(self.chainB_sc)

        self.num_sc = len(self.chainA_sc)

        self.patches_A_dict,self.patches_B_dict = self.assign_patches(self.protein)
        self.assign_biochem(self.protein)
    

    def assign_patches(self,protein):

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

        self.map_res_to_patchA = {}
        self.map_res_to_patchB= {}

        self.patch_labels ={'d1':[], 'd2':[], 'd3':[], 'M':[], 'NP-1':[], 'NP-2':[], 'NNP-1':[], 'NNP-2':[], 'tail':[], 'NNP-3':[], 'NNP-4':[]}
        self.patch_labelsB ={'d1':[], 'd2':[], 'd3':[], 'M':[], 'NP-1':[], 'NP-2':[], 'NNP-1':[], 'NNP-2':[], 'tail':[], 'NNP-3':[], 'NNP-4':[]}

        
        count=0
        for i in range(len(self.chainA_sc)):
            res = self.chainA_sc[i].resid

            if res in d1_A.resids:
                patches_A_dict[i]=1
                self.map_res_to_patchA[res]='d1'
                self.patch_labels['d1'].append(res)
            elif res in d2_A.resids:
                patches_A_dict[i]=2
                self.map_res_to_patchA[res]='d2'
                self.patch_labels['d2'].append(res)
            elif res in d3_A.resids:
                patches_A_dict[i]=3
                self.map_res_to_patchA[res]='d3'
                self.patch_labels['d3'].append(res)
            elif res in memb_A.resids:
                patches_A_dict[i]=4
                self.map_res_to_patchA[res]='M'
                self.patch_labels['M'].append(res)
            elif res in NN9_A.resids:
                patches_A_dict[i]=5
                self.map_res_to_patchA[res]='NP-1'
                self.patch_labels['NP-1'].append(res)
            elif res in patch2_A.resids:
                patches_A_dict[i]=6
                self.map_res_to_patchA[res]='NP-2'
                self.patch_labels['NP-2'].append(res)
            elif res in patch1_A.resids:
                patches_A_dict[i]=7
                self.map_res_to_patchA[res]='NNP-1'
                self.patch_labels['NNP-1'].append(res)
            elif res in patch3_A.resids:
                patches_A_dict[i]=8
                self.map_res_to_patchA[res]='NNP-2'
                self.patch_labels['NNP-2'].append(res)
            elif res in tail_A.resids:
                patches_A_dict[i]=9
                self.map_res_to_patchA[res]='tail'
                self.patch_labels['tail'].append(res)
            elif res in patch11_A.resids:
                patches_A_dict[i]=10
                self.map_res_to_patchA[res]='NNP-3'
                self.patch_labels['NNP-3'].append(res)
            elif res in patch13_A.resids:
                patches_A_dict[i]=11
                self.map_res_to_patchA[res]='NNP-4'
                self.patch_labels['NNP-4'].append(res)

        for i in range(len(self.chainB_sc)):
            res = self.chainB_sc[i].resid
            if res in d1_B.resids:
                patches_B_dict[i]=1
                self.map_res_to_patchB[res]='d1'
                self.patch_labelsB['d1'].append(res)
            elif res in d2_B.resids:
                patches_B_dict[i]=2
                self.map_res_to_patchB[res]='d2'
                self.patch_labelsB['d2'].append(res)
            elif res in d3_B.resids:
                patches_B_dict[i]=3
                self.map_res_to_patchB[res]='d3'
                self.patch_labelsB['d3'].append(res)
            elif res in memb_B.resids:
                patches_B_dict[i]=4
                self.map_res_to_patchB[res]='M'
                self.patch_labelsB['M'].append(res)
            elif res in NN11_B.resids:
                patches_B_dict[i]=5
                self.map_res_to_patchB[res]='NP-1'
                self.patch_labelsB['NP-1'].append(res)
            elif res in patch1_B.resids:
                patches_B_dict[i]=6
                self.map_res_to_patchB[res]='NP-2'
                self.patch_labelsB['NP-2'].append(res)
            elif res in patch2_B.resids:
                patches_B_dict[i]=7
                self.map_res_to_patchB[res]='NNP-1'
                self.patch_labelsB['NNP-1'].append(res)
            elif res in patch3_B.resids:
                patches_B_dict[i]=8
                self.map_res_to_patchB[res]='NNP-2'
                self.patch_labelsB['NNP-2'].append(res)
            elif res in tail_B.resids:
                patches_B_dict[i]=9
                self.map_res_to_patchB[res]='tail'
                self.patch_labelsB['tail'].append(res)
            elif res in patch9_B.resids:
                patches_B_dict[i]=10
                self.map_res_to_patchB[res]='NNP-3'
                self.patch_labelsB['NNP-3'].append(res)
            elif res in patch13_B.resids:
                patches_B_dict[i]=11
                self.map_res_to_patchB[res]='NNP-4'
                self.patch_labelsB['NNP-4'].append(res)

        # print(d1_A)
        # print(patch1_B)
        # print(self.patch_labels)
        # print(self.patch_labelsB)


        return(patches_A_dict,patches_B_dict)


    def assign_biochem(self,protein):

        self.biochem_dict={'Non-Polar':[], '+ Charge':[], 'Polar':[],'- Charge':[]}
        print(len(self.chainA_sc))
        for i in range(len(self.chainA_sc)):
            resname = self.chainA_sc[i].resname
            if resname == 'ARG' or resname == 'LYS' or resname == 'HIS':
                self.biochem_dict['+ Charge'].append(i)
            elif resname == 'ASP' or resname == 'GLU':
                self.biochem_dict['- Charge'].append(i)
            elif resname == 'SER' or resname == 'THR' or resname == 'TYR' or resname == 'ASN' or resname == 'GLN':
                self.biochem_dict['Polar'].append(i)
            else:
                self.biochem_dict['Non-Polar'].append(i)



if __name__ == '__main__':
    tpr_file = sys.argv[1]
    traj_file = sys.argv[2]
    PatchesClass(tpr_file,traj_file)
