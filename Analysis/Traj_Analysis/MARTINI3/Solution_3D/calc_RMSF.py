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
from MDAnalysis.core.universe import Merge


#Reading in all the inputs
struct_file = sys.argv[1]
traj_file = sys.argv[2]
outfile_pref = sys.argv[3]

#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()

#Selecting only protein beads
protein = u.select_atoms("byres name BB")
protein_rms = u.select_atoms("byres name BB")
protein_rmsf = u.select_atoms("byres name BB")
print("Total mass: ",protein.total_mass())
# sys.exit()
#print(protein)

#Selecting protein chains
pro_A = protein.select_atoms("bynum 1-469")
pro_B = protein.select_atoms("bynum 470-938")
pro_B_mic = protein_rms.select_atoms("bynum 470-938")
atom_len = len(pro_A)
pro_A_bb = pro_A.select_atoms("name BB")
pro_B_bb = pro_B.select_atoms("name BB")

n_frames = u.trajectory.n_frames
ts = u.trajectory.ts
print("Timestep:",ts.dt)
global box_dims
box_dims = ts.dimensions
print("Box dimension: ", box_dims)
print("Number of frames: ", n_frames)

"""
    --------------------------------------- RMSF -----------------------------------
    """

avg_struct = align.AverageStructure(u,u,select='bynum 1-469',ref_frame=0).run()
ref_struct = avg_struct.universe 
ref_group = ref_struct.atoms
# ref_group.write(outfile_pref+"_ChainA_average.pdb")
aligner_A = align.AlignTraj(u, ref_struct,select='bynum 1-469 and name BB',in_memory=True).run()

rmsf_chainA_sel = u.select_atoms('bynum 1-469 and name BB')
# rmsf_chainB_sel = u.select_atoms('bynum 470-938 and name BB')
rmsf_chainA = rms.RMSF(rmsf_chainA_sel).run()
# rmsf_chainBtoA = rms.RMSF(rmsf_chainB_sel).run()


"""
The calculation of RMSF for chainB is a little bit tricky. Ideally we would align only chain B
atoms to a reference frame first and then average the coordinates to get the Average Structure of chain B.

For example get the avg_structB = align.AverageStructure(u,u,select='bynum 470-938',ref_frame=0).run()
Then get the universe from this - ref_structB = avg_structB.universe
And then use this Universe for the Alignment to calculate RMSF. However there occurs a problem during this step
as the number of atoms in the reference structure(ref_structB) is only 469 (corresponding to chain B) which is not the same as the one in the trajectory.
Now during this alignment step we need to specify the atom selection to which the alignment is to be done. Here is where the problem lies - 
In the main trajectory atoms of chainB correspond to 470-938, but in the reference structure they correspond to 1-469. So it is not possible to align based upon this selection. 
The alignemnt can happend with the atom selection 'bynum 1-469 and name BB', but the RMSF calculation could be wrong since it is possible the AlignTraj function
is aligning the atoms of chain A to the reference structure of chain B!!!! ( I am not sure about this, but this is a possibility)

So to avoid this, I am calculating the average structure of chain B using the same selection as in the trajectory. This means it tries it best to align both chainA and chainB
to a ref frame and then average their coordinates. While this is not the ideal way to do it, the other method would require filtering teh trajectory to output only chain B atoms (from gromacs)
which could be a bit cumbersome. Also this method should work for the RMSF calculation. This is the following code: 

avg_structB = align.AverageStructure(u_copy,u_copy,select='bynum 1-938',ref_frame=0).run()
ref_structB = avg_structB.universe
ref_broupB = ref_structB.atoms
ref_broupB.write(outfile_pref+"_ChainB_average.pdb")
print(ref_structB.atoms)
aligner_B = align.AlignTraj(u_copy, ref_structB,select='bynum 470-938 and name BB',in_memory=True).run() 


But! But!. For calculating the RMSD w.r.t to an average structure, I think I need to output average structures specifically aligned to each chain. SO for that the code for chainB changes to - 
avg_structB = align.AverageStructure(u_copy,u_copy,select='bynum 470-938',ref_frame=0).run()
ref_structB = avg_structB.universe
ref_broupB = ref_structB.atoms
ref_broupB.write(outfile_pref+"_ChainB_average.pdb")

No point in calculating RMSF since To calculate the RMSF for chainB, the following command works - 
aligner_B = align.AlignTraj(u_copy, ref_structB,select='bynum 1-469 and name BB',in_memory=True).run()    #The ref structure for B only contains 469 atoms of chain B. 
But AGAIN this means that the atoms 1-469 in u_copy are aligned with atoms 1-469 in ref_structB. 
Unless this function has someway of validating atoms between the ref struct and teh traj (which I don't think is likely), it is definetly calculating the wrong RMSF

"""

avg_structB = align.AverageStructure(u_copy,u_copy,select='bynum 470-938',ref_frame=0).run()
ref_structB = avg_structB.universe
ref_broupB = ref_structB.atoms
# ref_broupB.write(outfile_pref+"_ChainB_average.pdb")
print(ref_structB.atoms)


merge_ref = Merge(ref_struct.universe.atoms,ref_structB.universe.atoms)
merge_ref.atoms.write(outfile_pref+"Averaged_Struct_chainAB.pdb")
# aligner_B = align.AlignTraj(u_copy, ref_structB,select='bynum 1-469 and name BB',in_memory=True).run()    #The ref structure for B only contains 469 atoms of chain B. 


# rmsf_chainB_sel = u_copy.select_atoms('bynum 470-938 and name BB')
# rmsf_chainB = rms.RMSF(rmsf_chainB_sel).run()
# rmsf_chainAtoB = rms.RMSF(rmsf_chainA_sel).run()

# outfileA = outfile_pref + "_RMSF_chainA.dat"
# with open(outfileA, 'a') as fl3:
#     fl3.write("#Residue\tRMSF\n")
#     for i in range(len(rmsf_chainA_sel.atoms)):
#         fl3.write("%d\t%.3f" %(i,rmsf_chainA.rmsf[i]))
#         fl3.write("\n")

        

#     fl3.write("\n")

# outfileB = outfile_pref + "_RMSF_chainB.dat"
# with open(outfileB, 'a') as fl3:
#     fl3.write("#Residue\tRMSF\n")
    
#     for i in range(len(rmsf_chainB_sel.atoms)):
#         fl3.write("%d\t%.3f" %(i,rmsf_chainB.rmsf[i]))
#         fl3.write("\n")

#     fl3.write("\n")