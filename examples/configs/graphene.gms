# Graphene cell
#include unitcells/graphene.gms

# Replicate the primitive cell
getin nx 5
getin ny 6
getin nz 1
+< reply $a$ 0   0   $nx$
+< reply 0  $b$  0   $ny$
+< reply 0   0  $c$  $nz$

# Modify the lattice metrics.
getin r0 1.421437
set expand $r0$ $r0$ 1
box expand $r0$ $r0$ 1

# Set the element name
set element C

# Coordinates offset
set cm_pos 0 0 0 
          
# Complete borders with hydrogen atoms
getin r1 {1.2*$r0$}
getin r2 {1.3*$r0$}
- bo_range -1 0.9 1.0 $r1$ $r2$
set sp 2
+< fillh
          
# Save the coordinates
out posxyz graphene_ribon.xyz
                              

 
# Graphene cell
#include unitcells/graphene.gms

# Replicate the primitive cell
getin nx 30
getin ny 30
getin nz 1
+< reply $a$ 0   0   $nx$
+< reply 0  $b$  0   $ny$
+< reply 0   0  $c$  $nz$

# Modify the lattice metrics.
getin r0 1.421437
set expand $r0$ $r0$ 1
box expand $r0$ $r0$ 1

# Set the element name
set element C

# Coordinates offset
set cm_pos 0 0 0 
                       
# Cutting and rotating to form a triangle
set cm_pos 0 0 0 
set rotate z  120
>> xrange -80  15
set rotate z  120
>> xrange -80  15
set rotate z  120
>> xrange -80  15
   
# Complete borders with hydrogen atoms
getin r1 {1.2*$r0$}
getin r2 {1.3*$r0$}
- bo_range -1 0.9 1.0 $r1$ $r2$
set sp 2
+< fillh
          
# Save the coordinates
out posxyz graphene_triangle.xyz
         
               
