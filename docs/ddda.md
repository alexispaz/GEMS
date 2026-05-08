# Gráfico de Flyvbjerg-Petersen.

La varianza de 

Tomando un numero finito de muestras $N$ sobre una cantidad estadistica $X$, la
media define un estimador para el valor de expecation $\langle X \rangle$:

$\bar{X}=\sum_i X_i/N$
 
Se suele reportar el resultado de una medicion sobre $X$ como $(\bar{x} +-
E_s)$, siendo $E_s$ el error estandar:

$E_s=t\sigma/\sqrt(N)$
 
Donde $t$ es el factor de confianza (i.e. 1.92 para el 95% de confianza)
para el cual se considera el error y $\sigma$ es la desviación estandar de $X$
dada por:

$\sigma^2=\bar{X^2}-\bar{X}^2$

Usando $\bar{X}$ en la expresion anterior, se puede define un estimador para
$\sigma$ para las $N$ muestras. Al estimar $\sigma$, notesé que $E_s$ tambien
adquiere su propio error estandar [^1]:

$\sigma_{E_s}=E_s \sqrt{2/(N-1)}$

Esto supone que las $N$ muestras no está correlacionadas. Si se toman $N$
muestras que si están correlacionadas, se puede tomar la media por bloques de
$m>1$ muestras y considerar las resultantes como muestras. Estas $N/m$ muestras
tienden a estar decorrelacionadas a medida que $m$ se acerca a $N$.  Un grafico
típico de del estimador de $\sigma$  con su error en funcion del $log_2(m)$
sería [^1]:

    0.0013 ++----------------------------------------------------------------+
           |                                                    ***          |
    0.0012 ++                                               ***  *      ***  |
    0.0011 ++                                                *   *       *   |
           |                                       ***       *   *  ***  *   |
     0.001 ++   Es=sqrt(var/N-1)           *** ***  A  ***   A   *   *   *   |
           |                           ***  A   A   *   A    *   A   *   *   |
    0.0009 ++                          *A* *** *** ***  *    *   *   *   *   |
           |                       *A*                 ***  ***  *   *   *   |
    0.0008 ++                                                    *   A   A   |
    0.0007 ++                  *A*                              ***  *   *   |
           |                           barras de error:              *   *   |
    0.0006 ++                            Es/(sqrt(2(N-1))            *   *   |
           |               *A*                                       *   *   |
    0.0005 ++                                                       ***  *   |
           |          *A*                                                *   |
    0.0004 ++                                                           ***  |
    0.0003 ++     *A*                                                        |
           +  *A*  +        +       +       +       +        +       +       +
    0.0002*A*------+--------+-------+-------+-------+--------+-------+-------+
           0       2        4       6       8       10       12      14      16
 
 
Un plato en este grafico indica que se alcanzo un numero de muestras suficiente
decorrelacionadas para estimar $\sigma$ correctamente. El algorithmo DDDA
propuesto por Kent et. al.[^2] permite calcular $\sigma^2$ este gráfico
"on-the-fly" en una simualción. Se encuentra implementado en el modulo
`gems_ddda`. 

Un criterio para definir "on-the-fly" el plato de este gráfico fue propuesto en
[^3] y se encuentra implementado en la subrrutina `decorrelation_variance` de
este módulo. La idea es partir del bloque con $m=N$ e ir intersectando el
intervalo de error hacia los bloques mas pequeños. Cuando la interseccion es
vacia se consdiera el número de puntos/barras de error que se intersectaron en
comparacion con el número de puntos totales del gráfico.

# Dynamic Distributable Decorrelation Algorithm (DDDA)

Este algorithmo fue propuesto por Kent et. al.[^3] permite calcular
$\sigma^2$ "on-the-fly".

DDDA requiere O(N+mlog_2(N)) operaciones para calcular la varianza m veces
durante una simulacion de N paso Si se programa parallelo requiere O(log_2N)
de datos a comunicar y el bond procedure. 

A medida que el dato es generado durante el calculo, es sumado a un objeto de
la clase decorrleation. Hay que tener en cuenta que el dato se pierde, en el
sentido de que no se puede retirar despues (no se puede hacer un buffer que
vaya avanzando)

Calculos en paralelo: Si la idea es que funcione de manera serial, solo un
objeto decorrelation basta. Si uno quiere que funcione de manera paralelo, cada
procesador, que genera sus datos, debe tener su propio objeto de decorrelacion
y luego se los suma a travez del proceso decorrelation_addition para obtener la
desviacion estandar de propiedad calculada en paralelo.

---

[^1]: Flyvbjerg, H. and Petersen, H. G. (1989), Error estimates on averages of
correlated data. J. Chem. Phys. 91, 461-466. https://doi.org/10.1063/1.457480

[^2]: Paz, S. A. and Leiva E. P. M. (2015) Time Recovery for a Complex Process
Using Accelerated Dynamics. J. Chem. Theory Comp. 11 (4), 1725-1734
http://doi.org/10.1021/ct5009729
 
[^3]: Kent, D.R., IV, Muller, R.P., Anderson, A.G., Goddard, W.A., III and
Feldmann, M.T. (2007), Efficient algorithm for “on-the-fly” error analysis of
local or distributed serially correlated data. J. Comput. Chem., 28
2309-2316. https://doi.org/10.1002/jcc.20746
 
