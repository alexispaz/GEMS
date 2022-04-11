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
        
# Define a sistem
box size 5.78 5.78 5.78

>< atom 1 1 1
+< atom 2.2 1 1
# m0=1uma
sys add H  

> sys
set pbc F F F 
set mass 3.
group 1 add

# Get variables
> atom 1
getin x1 cm_pos
> atom 2
getin x2 cm_pos
       
# Testing interaction commands
interact 1 pair lj 1. 0.890898718140339 2

# Testing math
tau_ps:=0.01018
time step {0.002*$tau_ps$}   (integration timestep [ps])

# Testing evolve commands  
> sys
evolve ermak 100 1
evolve :1 v_verlet
bloque load_each 3 2
dinamica 6
  
# Print properties
> atom 1
print cm_pos
> atom 2
print cm_pos
  
# Testing set/get
> atom 1
set cm_pos $x1[1] $x1[2] $x1[3]
print cm_pos
 
