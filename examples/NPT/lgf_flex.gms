(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

# This LJ parameters was used in:  
# Grønbech-Jensen, N., & Farago, O. (2014). Constant pressure and
# temperature discrete-time Langevin molecular dynamics. The Journal of
# Chemical Physics, 141(19), 194108. http://doi.org/10.1063/1.4901303
                     
dimension 3
 
# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
prng seed 123456
 
#include ../configs/unitcells/fcc.gms
+< reply 1 0 0   6
+< reply 0 1 0   6
+< reply 0 0 1   6
d:=1.
set expand {$d$*sqrt(2)} {$d$*sqrt(2)} {$d$*sqrt(2)}
box size 6 6 6
box expand {$d$*sqrt(2)} {$d$*sqrt(2)} {$d$*sqrt(2)}

# m0=1uma
set element H


> sys
  set pbc T T T 
  group 1 add
box mic F

((((((((((( caja )))))))))))

box move (make the total velocity zero)
> sys
set move 0.2 0.2 0.2

# Epsilon=1eV=1E0
# Rmin=1A=1r0
sigma:={2**(-1./6.)}
interact 1 pair slj 1. $sigma$ 1.108683
               
# About time unit
#   1eV=9648.61 uma*A**2/ps**2
#   1E0=9648.61 m0*r0**2/ps**2
#   1ps**2=  9648.61 m0*r0**2/E0
#   1ps=sqrt(9648.61) tau 
#   1ps=98.22 tau
tau_ps:=0.01018
       
time step {0.01*$tau_ps$}   (integration timestep [ps])
                                         

((((((((((( archivos de salida  )))))))))) 

outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1 temp 1
outfile :f1 each 10
# This should be in 1
# plot [][:] 'Energy.gms.dat' u ($0*100*0.002):($4*kB_eV) w l 
           
outfile :f2 name Press.$jobname$.dat
outfile :f2 cols pressure 1
outfile :f2 each 10
               
outfile :f3 name Pos.$jobname$.xyz
outfile :f3 pos 1
outfile :f3 each 10
         
outfile :f4 name Caja.$jobname$.dat
outfile :f4 cols box
outfile :f4 each 10
            
outfile :f5 name Viri.$jobname$.dat
outfile :f5 cols virial 1
outfile :f5 each 10
           

> sys

# T=1E0/kB=eV/kB=11604.45 K
kt_k:=11604.45

# set tempgdist {0.3*$kt_k$}
# add v_verlet 

# add ermak {0.3*$kt_k$} {1/$tau_ps$} 
# out state

# Note that alpha is gama*m

# P=1eV/A^3=1E0/r0
e0r0_atm:=1581225.25240563
evolve piston_lgf x {0.3*$kt_k$} {1/$tau_ps$} {0.1*$e0r0_atm$} {1/$tau_ps$} 0.1
evolve piston_lgf y {0.3*$kt_k$} {1/$tau_ps$} {0.1*$e0r0_atm$} {1/$tau_ps$} 0.1
evolve piston_lgf z {0.3*$kt_k$} {1/$tau_ps$} {0.1*$e0r0_atm$} {1/$tau_ps$} 0.1


dinamica 100
