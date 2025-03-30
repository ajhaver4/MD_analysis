#TCL Scripting

mol new solution_protein.pdb
mol new lsp_welltemp_244us.gro

mol selection "chain A"
mol addrep 0
mol modstyle 1 0 NewCartoon
mol modcolor 1 0 ColorID 19
mol modmaterial 1 0 AOShiny

mol selection "chain B"
mol addrep 0
mol modstyle 2 0 NewCartoon
mol modcolor 2 0 ColorID 31
mol modmaterial 2 0 AOShiny

mol selection "chain A"
mol addrep 0
mol modstyle 3 0 QuickSurf 1.5 0.8 0.8 1.0
mol modcolor 3 0 ColorID 19
mol modmaterial 3 0 Glass1

mol selection "chain B"
mol addrep 0
mol modstyle 4 0 QuickSurf 1.5 0.8 0.8 1.0
mol modcolor 4 0 ColorID 31
mol modmaterial 4 0 Glass1

#Ions display
mol selection "name NA or name CL"
mol addrep 1
mol showrep 1 0 off
mol modstyle 1 1 Beads 2.5 5
mol modmaterial 1 1 Glass3
mol modcolor 1 1 ColorID 0

#Turn off 'all' rep
mol showrep 0 0 off

#Display Settings
color Display Background white
axes location off
display projection orthographic
display depthcue off
display antialias on
display aoambient 0.4
display aodirect 0.6
light 0 on
light 1 off
light 2 on
light 2 pos {5 2 -1}
light 3 on


#Translation and rotations





#



