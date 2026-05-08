# Input script for GEMS

# Testing blocks
bloque repeat 12
  print {2*$i$}
fin

bloque save 1
   print hola
fin
    
bloque save 2
   print dm_steps
fin
    
bloque load 1
        
# Box
box size 5.78 5.78 5.78

# Play with particles that will not be used
+ atom 2.2 1 10
> all
^+ reply 0 0 1 10
^- zrange 15 20
^~ zrange 13 20
group 2 add

# Particles that will be used
>+ atom 1 1 1
^+ atom 2.2 1 1
set element H

> element H
set pbc F F F 
set mass 3.
group 1 add

# Save variables for later use
^> atom 1
getin x1 cm_pos
> group 1
^> atom 2
getin x2 cm_pos

# Testing interaction commands
interact 1 pair lj 1. 0.890898718140339 2

# Testing math
tau_ps:=0.01018
time step {0.002*$tau_ps$}   (integration timestep [ps])

# Testing evolve commands  
> all
evolve ermak 100 1
evolve :1 v_verlet
bloque load_each 3 2
dinamica 6
  
# Print properties
> atom 7
print cm_pos
> atom 8
print cm_pos
  
# Testing set/get
> group 1
^> atom 1
set cm_pos $x1[1] $x1[2] $x1[3]
print cm_pos
 
# Final deletion of particles
- group 2

