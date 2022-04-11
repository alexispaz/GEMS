(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3

box size    13.1399  12.64396  20.0000 

>+ read Pos.md.xyz

> all
group 1 add

> all
set pbc t t t

interact 1 pair lj 0.1 3.890898718140339 10.

outfile :f3 name Energy.$jobname$.dat
outfile :f3 cols epot 1
outfile :f3 each 1
    
> all

dinamica_from_xyz 6 Pos.md.xyz
