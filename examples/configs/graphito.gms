# Input script for GeMS

# This script creates a graphite crystal. General procedures for creation of
# orthorombic system can be found in the orthorhombic_crystal.gms example

# Graphene cell
#include unitcells/graphene.gms

# Replicate the primitive cell
nx:=3
ny:=6
nz:=3
^+ reply $a$ 0   0   $nx$
^+ reply 0  $b$  0   $ny$
^+ reply 0   0  $c$  $nz$

# Modify the lattice metrics.
set expand 1.421437 1.421437 3.35
box expand 1.421437 1.421437 3.35  

# Coordinates offset
set minpos x 0.5
set minpos y 0.5
set minpos z 0.5

# Set the element name
set element C

# Save the coordinates
out posxyz $jobname$.xyz

