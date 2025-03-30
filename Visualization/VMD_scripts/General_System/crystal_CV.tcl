
#Load LSP dimer
mol new final_chainAB.pdb

#Display Settings
color Display Background white
axes location off
display projection orthographic
display depthcue off
display antialias on
display aoambient 0.4
display aodirect 0.6
light 0 on

#Turn off other rep 
mol showrep 0 0 off


#Change representations
mol selection "chain A"
mol addrep 0
mol modstyle 1 0 NewCartoon
mol modcolor 1 0 ColorID 0
mol modmaterial 1 0 AOShiny

mol selection "chain B"
mol addrep 0
mol modstyle 2 0 NewCartoon
mol modcolor 2 0 ColorID 28
mol modmaterial 2 0 AOShiny


mol selection "resid 57 to 64 and chain A"
mol addrep 0
mol modstyle 3 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 3 0 ColorID 14

mol selection "resid 78 to 92 and chain B"
mol addrep 0
mol modstyle 4 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 4 0 ColorID 14

mol selection "resid 224 to 229 and chain A"
mol addrep 0
mol modstyle 5 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 5 0 ColorID 19

mol selection "resid 189 to 196 and chain B"
mol addrep 0
mol modstyle 6 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 6 0 ColorID 19
