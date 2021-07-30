
dimension 3

# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
prng seed 123456

>< read conf.xyz
sys add

> sys
set pbc F F F 
group 1 add
> element Co
set mass 40.
group 2 add 
> element Au  
group 3 add
set mass 40.



box move (make the total velocity zero)

interact 2 pair lj  0.10364 3.4 200
interact 3 pair lj  0.06218 3.4 200
interact 2 < 3 pair lj 0.08291 3.4 200
interact 3 < 2 pair lj 0.08291 3.4 200

time step 0.005   (integration timestep [ps])
                
                        
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
outfile :f1 each 1e6

outfile :f2 name Pos.$jobname$.xyz
outfile :f2 pos 1
outfile :f2 each 1e6
        

> sys
out state 

> sys
set tempgdist 500
evolve ermak 500 2     


setwallauxcore  0.0 60.0
setwallauxshell 0.0 60.0
setwall1d 0.0 2.0 100 17.5
setwtmd1d 0.1 0.01 9000 0.1 1e5 13 500
wtdcm 1e3

outfile :f3 name E_libre.dat
outfile :f3 free_en_1d 1
outfile :f3 each 1e6
 
out state 
