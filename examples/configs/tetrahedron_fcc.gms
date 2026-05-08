# Input file for GeMS

# Create a tetrahedron from cutting a fcc cube. 
# All the resulting faces are (111).

# Create the cube
#include unitcells/fcc.gms
^+ reply   1 0 0     4
^+ reply   0 1 0     4
^+ reply   0 0 1     4

out posxyz asd.xyz

# Perform the cutting operation using planes
^> below   3.5  3    3     1  2.5  0.5   1  0.5 2.5
^> below   3.5  3    3     1  2.5  0.5   3  0.5 0.5 
^> below   3.5  3    3     1  0.5  2.5   3  0.5 0.5
^> above   1   2.5  0.5    1  0.5  2.5   3  0.5 0.5

# Rotations to align a face on the xy-plane.
set rotate z -45    3 0.5 0.5
set rotate y -54.74 3 0.5 0.5
set move 0 0 -0.5

# The tip is at 2.183326 1.914214 2.309339
out posxyz $jobname$.xyz  

