(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3

>< read ../configs/ref/gold_cube.xyz
sys add Au

# box size 12.19500   8.13000  12.19500 
box size 12.19500   15  12.19500 

>< atom {$boxx$*0.2} 9.13000 {$boxz$*0.2}
sys add Pt

> sys 
group 1 add
set pbc f f f
set move 0.5 0.5 0.5

> element Pt
group 2 add

interact 1 under tb read test.prm
 
# Outfiles
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols time energy 1
outfile :f1 each 1
    
outfile :f2 name Pos.$jobname$.xyz
outfile :f2 Poscr 1
outfile :f2 each 1


# Image number and system dimension
> sys
neb_nimg 30

# Initial configuration
> sys
lbfgs
neb_img 1

# Middle configuration
> group 2
set move 2 0 2
> sys
lbfgs
neb_img 15
          
# Final configuration
> group 2
set move 2 0 -2
> sys
lbfgs
neb_img 30
          
 
# 10 iterations with a spring constant of 5 
neb_run 20 5
                            
