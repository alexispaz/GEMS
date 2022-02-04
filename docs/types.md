# Atoms and groups

There are two important classes in GEMS

atom
: An **atom** is an atom. It has some basic components used by a large number
of algorithms, like the position, mass, chemical specie or velocity.

group
: A **group** is a group of **atoms**. In general, any algorithm is applied to
a certain group. The definition of a group allows to mix different algorithms.

Beyond the basic components of an **atom**, different algorithms may require
other variables. For instance, a simple Lennard-Jones potential may require to
define two extra-parameters (e.g. epsilon and sigma) or more. Another example,
a 5th-order Predictor-Corrector integration algorithm will need to track more
derivatives of the atom coordinates than a simple velocity Verlet algorithm.
To store all this variables in a simple **atom** class will be very
inefficient, i.e. the class size will grow with the number of algorithm in the
code. On the other hand, if we create all the **group** class procedure using a
parent **atom** and then use different **atom** childs for each algorithm,
`select type` blocks will be required during the algorithm iteration.
 
One way is to create a parent **group** class with different **group** childs for each
algorithm that holds the extra atomic variables. If we attempt to build
dynamic memory structures for this variables using their intrinsic types
(i.e. structures for `integer`, `real`, etc.), it might be unpractical to
associate each variable to its **atom**s. For instance, considering just an
example with two extra variables, the group iteration will become something
like:

	la=>group%alist
	ls=>group%sigma
	le=>group%epsilon
	do i=1,group%nat
	  la=>la%next
	  ls=>ls%next
	  le=>le%next
	enddo

And if we group all them in user defined types, like this:

	la=>group%alist
	lp=>group%properties
	do i=1,group%nat
	  la=>la%next
	  lp=>lp%next
	  sigma=lp%s
	  epsilon=lp%e
	enddo

But we need to code the `group%properties` dynamic structure procedures for
each user defined type. We can attempt to introduce a
class pointer from the `alist` to this kind of property structure, but then
`select type` blocks appears again.

So, I think plain arrays in extended **group** types may be the best. Although
these arrays should be reallocated when atoms are added to the group, an
arbitrary grow step can be used to avoid a low performance. In a similar way,
memory rearrangements due to atom deletions can be fix made in a
step fashion, for instance, when more than 100 atoms are deleted. However, to
associate each **atom** with its corresponding array element, we needed to
create and maintain and index. Still, for those simple algorithms
that do not require further **atom** parameters, using a basic group
object without index may be better. Therefore, I keep two kind of groups:

- A **group** class that use a dynamic linked list of atoms.

- An **igroup** that keeps an index of atoms and can be extended to include
arrays with new atom properties.
 
To cyle a **group** we can use:

    la => g%alist
    do ii = 1,g%ref%nat
      la => la%next

I make the **igroup** a child of the **group**, retaining also the linked list
to re use the **group** procedures. So we can cycle an **igroup** in the same
way that we cycle a **group**. From the **igroup** structure, we can browse
the index of atoms using an array of pointers. On the other hand, an
array of integers, probably small, can be included in each **atom** to access
its ids, i.e. its index positions in each **igroup** it belongs. We can cycle
the **igroup** and access to the index in two ways. First:

    la => g%alist
    do ii = 1,g%nat
      la => la%next
      a => la%o
      i = a%gid(g)

where `gid` is a function that find the in the internal index of the atom.
Second:
       
    do i = 1,g%amax
      a => g%a(i)%o
      if(.not.associated(a)) cycle

Where `g%a` is the array of pointers and the conditional filter out already
detached atoms. `g%amax` is always grather or equal to `g%nat` and its
difference can be used to decide when to update the structure, or even to
reduce the arrays sice in case memory is important.

There are two particular extensions of igroups that are worth to discuss
more.

- The **cgroup** includes components and procedures for sorting atoms in cells. 

- The **ngroup** includes a neighbor list and a searching procedure. It has also two **group**
components inside: `ngroup%b` holds the potential atom neighbors of the atoms in the
`ngroup%ref`. `ngroup%b` is a **cgroup** while `ngroup%ref` a simple **group**. 

