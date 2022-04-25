(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    Estableciendo la  configuracion               )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3
box size 50 50 50

# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
prng seed 123456
 
>+ fill 100 2.6
set element He

> all
set pbc T T T 
group 1 add
 


((((((((((( caja )))))))))))


box move (make the total velocity zero)
time step 0.002d0  (integration timestep [ps])

((((((((((( interacciones  )))))))))) 

# ==============================================================
# Taked from par_all22_prot.inp CHARMM parameter file:
# 
# NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
# cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5 
#                 !adm jr., 5/08/91, suggested cutoff scheme
# 
# V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
# epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
# Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
# HE     0.000000  -0.021270     1.4800  ! helium, experimental pot. energy surface, adm jr., 12/95
# ==============================================================

cal_ev := 2.61144768e19
na := 6.0221417930e23
kcalmol_ev := {$cal_ev/$na*1000.}
interact 1 pair plj {0.021270*$kcalmol_ev} 1.4800

((((((((((( archivos de salida  )))))))))) 

outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1 temp 1
outfile :f1 each 100

# outfile :f2 name Pos.$jobname$.xyz
# outfile :f2 pos 1
# outfile :f2 each 100
 
> all
set tempgdist 300  
evolve :pepe ermak 300 1
evolve gcmc  {$boxz*4./5} {$boxz*5./5} 0.005 1. 10 300.
out state
dinamica 1000

