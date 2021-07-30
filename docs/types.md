# Atoms and groups

There are two important classes in GEMS: the **atom** and the **group**.

- An **atom** is an atom. It has some basic components used by a large number
of algorithms, like the position, mass, chemical specie or velocity.

- A **group** is a group of **atoms**. In general, any algorithm is applied to
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

And if we group all them in user-types, like this:

	la=>group%alist
	lp=>group%properties
	do i=1,group%nat
	  la=>la%next
	  lp=>lp%next
	  sigma=lp%s
	  epsilon=lp%e
	enddo

Then we need to code the `group%properties` dynamic structure procedures for
each child **group**. We can attempt to introduce a class pointer from the
`alist` to this kind of property structure, but then `select type` blocks
appears again.

So, I think plain arrays may be the best. If atoms are added or deleted from
a group, these arrays should be reallocated. To associate each **atom**
with its arrays, an internal **group** index is needed. Still, for those algorithms
that do not require further parameters than the basic components of the
**atom**, dynamic memory groups may be better. Therefore, I
keep two kind of groups:

- A **group** class that use a dynamic linked list of atoms.

- An **igroup** that keeps an index of atoms and can be extended to include
arrays with new atom properties.

I make the **igroup** a child of the **group**, retaining also the linked list
to re use the **group** procedures. The atom index can be maintained using an
array of pointers in the **igroup**. An array of integers can be included in
each **atom** to save the index associated to the **igroup** it belong. When an
**atom** pointer is attached to a **group**, this **group** is also added to
this array of ids.

There are to extensions of igroups that are worth to discuss more.

- The **cgroup** includes components and procedures for sorting atoms in cells. 

- The **ngroup** includes a neighbor list and a searching procedure. It has also two **group**
components inside: `ngroup%b` holds the potential atom neighbors of the atoms in the
`ngroup%ref`. `ngroup%b` is a **cgroup** while `ngroup%ref` a simple **group**. 

