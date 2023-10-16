# shells
Programs for computing the redistribution of lipids around a protein

# generate borders.bin file (put in your psf/pdb file as well as a time-ordered series of dcd files: 
hopBorders psf dcd

# process the borders file, the threshold 0.5 here will ignore small borders
processHop borders.bin 0.5 > shells.txt

# either a 1 or -1 for each lipid. here for example everything less than or equal to res 340 was lower:
grep " C1 " frame.pdb | awk '{ if( $5 <= 340 ) print " -1"; else print " 1"; }' | tr -d '\n' > leaflet.txt

# gather the lipid index from a frame of the pdb
outputIndex frame.pdb > lipid.index

# collect a histogram from shells.txt
histo lipid.index shells.txt --leaflet=leaflet.txt --upper

# histograms can be averaged with avhist
