# Input script for GeMS

# This script creates a tetragonal crystal bulk using an already defined
# primitive cell.
          
# Uncomment the desired primitive cell:
# 
# Face-centered diamond-cubic  
# include unitcells/fcdc.gms
# 
# Face-centered cubic
#include unitcells/fcc.gms
# 
# Graphene
# include unitcells/graphene.gms

# Optional: modify the lattice point metrics. Intended for crystals that have
# complex lattice points (e.g. molecules, nanoparticles, etc.).  
# set expand 1 1 1
         
# Optional: modify the lattice parameters. Included primitive cell have defined
# the variables a, b and c as the primitive cell size. Intended for crystals
# that have complex lattice points (e.g. molecules, nanoparticles, etc.).
# getin a 3
# getin b 3
# getin c 3

# Replicate the primitive cell
getin nx 3
getin ny 2
getin nz 3
+< reply $a$ 0   0   $nx$
+< reply 0  $b$  0   $ny$
+< reply 0   0  $c$  $nz$
               
# Update the box size
box size   $nx$ $ny$ $nz$
box expand $a$ $b$ $c$
              
# Optional: modify the lattice metrics. Useful to change neighbor lengths for
# simple crystals.
# 
# For graphite:
# set element C
# set expand 1.421437 1.421437 3.35
# box expand 1.421437 1.421437 3.35
# 
# # For c-Si:
# set element Si
# set expand 5.43071 5.43071 5.43071
# box expand 5.43071 5.43071 5.43071
# 
# For gold:
set element Au
set expand 4.065 4.065 4.065
box expand 4.065 4.065 4.065
             
# Save the configuration
out posxyz gold_cube.xyz
