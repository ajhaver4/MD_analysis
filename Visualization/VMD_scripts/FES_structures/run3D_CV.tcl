
#Load LSP dimer
mol new Alignedd1_Avg_3D.gro_atomistic.gro

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
mol selection "serial 1 to 3396"
mol addrep 0
mol modstyle 1 0 NewCartoon
mol modcolor 1 0 ColorID 0
mol modmaterial 1 0 AOShiny

mol selection "serial 3397 to 6792" 
mol addrep 0
mol modstyle 2 0 NewCartoon
mol modcolor 2 0 ColorID 28
mol modmaterial 2 0 AOShiny


mol selection "serial 92 to 195 or serial 218 to 231" 
mol addrep 0
mol modstyle 3 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 3 0 ColorID 14

mol selection "serial 3826 to 4064"
mol addrep 0
mol modstyle 4 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 4 0 ColorID 14

mol selection "serial 2787 to 2881"
mol addrep 0
mol modstyle 5 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 5 0 ColorID 19

mol selection "serial 5611 to 5759"
mol addrep 0
mol modstyle 6 0 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 6 0 ColorID 19
