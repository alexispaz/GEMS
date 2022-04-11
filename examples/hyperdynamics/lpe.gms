(((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))
(                        INPUT SYSTEM FOR GCMD PROGRAM                         )
(                  comentaries are between brackets or whit "# "               )
(                          caps can be lower or upper                          )
(                   caps must be upper for potential names                     )
(((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))

dimension 3
 
# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
prng seed 123456
 
>< read clst_a42.xyz

# Tomo una cara aleatoria
set cm_pos 0.0 0.0 0.0

set rotate x {$rnd*360}
set rotate y {$rnd*360}
set rotate z {$rnd*360}

set maxpos x 10
set minpos y 1
set minpos z 1
set element Au

>< read clst_Ih0013i.xyz
set minpos x 13
set minpos y 1
set minpos z 1
set element Co
 
box move

time step 0.005

> sys
group 1 add

> element Au
group 2 add

> element Co
group 3 add

interact 1 tb read ../../parameters/lucas_2017.prm
interact :hybrid_ea 1 bias compress_below 1.0 1.0
    

# > element Au 
# set move 100 100 100
# lbfgs

# > element Co 
# lbfgs

# > element Au 
# set move -100 -100 -100

(selecciono la salida)
> group 1


outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols time step epot 1 bias
outfile :f1 each 100

outfile :f2 name Pos.$jobname$.xyz
outfile :f2 pos 1
outfile :f2 each 100
 
outfile :f4 name FPp.$jobname$.dat
outfile :f4 hd_fpp
outfile :f4 each 100
outfile :f4 at hd
 
> sys
set tempgdist 300  
evolve :termo ermak 300 10.
evolve :termo voter .true.

ngrad:={55*3} ( numero de grados de libertad )
nblock:=7     ( numero de bloques            ) 
nsteps:=20000 ( pasos por bloque             ) 
nequil:=1000  ( pasos de equilibracion       ) 
ep:=0.3       ( eprime o w                   ) 
ap:=50.       ( aprime o B                   ) 
plato:=0.36   ( factor plato                 ) 
pbmax:=0.005  ( probabilidad bias de corte   ) 
temp:=300.0   ( temperatura                  ) 

hybrid_wb   $ngrad $nblock $nsteps $nequil $ep $ap $plato $pbmax $temp
# dinamica 100

