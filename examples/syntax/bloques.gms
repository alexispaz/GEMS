# Input script for GeMS

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
        
box size 5.78 5.78 5.78

>< atom 1 1 1
+< atom 2.2 1 1
# m0=1uma
sys add H  
 
> sys
set pbc F F F 
set mass 3.
group 1 add
      
interact 1 pair lj 1. 0.890898718140339 2

tau_ps:=0.01018
time step {0.002*$tau_ps$}   (integration timestep [ps])

> sys
evolve ermak 100 1
evolve :1 v_verlet
bloque load_each 3 2
dinamica 6
