
#Load LSP dimer
mol new Alignedd1_Avg_2D.gro_atomistic.gro

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
mol showrep 1 0 off


#Change representations
mol selection "serial 1 to 3396"
mol addrep 1
mol modstyle 1 1 NewCartoon
mol modcolor 1 1 ColorID 0
mol modmaterial 1 1 AOShiny

mol selection "serial 3397 to 6792" 
mol addrep 1
mol modstyle 2 1 NewCartoon
mol modcolor 2 1 ColorID 28
mol modmaterial 2 1 AOShiny


mol selection "serial 92 to 195 or serial 218 to 231" 
mol addrep 1
mol modstyle 3 1 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 3 1 ColorID 14

mol selection "serial 3826 to 4064"
mol addrep 1
mol modstyle 4 1 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 4 1 ColorID 14

mol selection "serial 2787 to 2881"
mol addrep 1
mol modstyle 5 1 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 5 1 ColorID 19

mol selection "serial 5611 to 5759"
mol addrep 1
mol modstyle 6 1 QuickSurf 1.0 0.8 0.8 1.0
mol modcolor 6 1 ColorID 19
