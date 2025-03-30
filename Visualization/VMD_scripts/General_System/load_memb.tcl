mol new rearrange_protein_onmembrane.pdb
mol new memb_only.gro

mol selection "serial 1 to 3396"
mol addrep 0
mol modstyle 1 0 NewCartoon
mol modcolor 1 0 ColorID 19
mol modmaterial 1 0 AOShiny

mol selection "serial 3397 to 6792"
mol addrep 0
mol modstyle 2 0 NewCartoon
mol modcolor 2 0 ColorID 31
mol modmaterial 2 0 AOShiny

mol selection "serial 1 to 3396"
mol addrep 0
mol modstyle 3 0 QuickSurf 1.5 0.8 0.8 1.0
mol modcolor 3 0 ColorID 19
mol modmaterial 3 0 Glass1

mol selection "serial 3397 to 6792"
mol addrep 0
mol modstyle 4 0 QuickSurf 1.5 0.8 0.8 1.0
mol modcolor 4 0 ColorID 31
mol modmaterial 4 0 Glass1


#Membrane
mol selection "name P4 or name P5 or name PO4 or name CNO or name NC3 or name C1 or name C2 or name C3"
mol addrep 1
mol modstyle 1 1 QuickSurf 1.5 0.6 1.2 1.0
mol modcolor 1 1 ColorID 14
mol modmaterial 1 1 Glass2

mol selection "name C1A or name C1B or name C2A or name C2B or name C3A or name C3B or name C4 or name C4A or name C4B or name D2A"
mol addrep 1
mol modstyle 2 1 QuickSurf 1.2 1.0 1.3 1.0
mol modcolor 2 1 ColorID 6
mol modmaterial 2 1 Glass2

mol selection "name GL1 or name GL2"
mol addrep 1
mol modstyle 3 1 QuickSurf 1.2 1.0 1.3 1.0
mol modcolor 3 1 ColorID 17
mol modmaterial 3 1 Transparent


#Turn off 'all' rep
mol showrep 0 0 off
mol showrep 0 1 off

#Display commandscolor Display Background white
axes location off
display projection orthographic
display depthcue off
display antialias on
display ambientocclusion on
display aoambient 0.4
display aodirect 0.6
light 0 on
light 1 off
light 3 on
light 3 pos {10 10 0.0}

