
#Load LSP dimer
mol new ./MARTINI2/lspAB_cry.gro_atomistic.gro

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
mol showrep 3 0 off


#Change representations
mol selection "serial 1 to 3396"
mol addrep 3
mol modstyle 1 3 NewCartoon
mol modcolor 1 3 ColorID 0
mol modmaterial 1 3 GlassBubble

mol selection "serial 3397 to 6792" 
mol addrep 3
mol modstyle 2 3 NewCartoon
mol modcolor 2 3 ColorID 28
mol modmaterial 2 3 GlassBubble


