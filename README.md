nkmax, nkith, nkithr, popldump
==============================

`nkmax` generates a population of *n* unsigned integers distributed among
*k* nodes and finds the maximum element.

`nkith` and `nkithr` both generate *n* unsigned integers distributed among 
*k* nodes and find the *i*th element.

`popldump` dumps the binary debug output file of `nkmax`, `nkith`, and `nkithr`
in a human readable format.

NB: `nkith` is not recommended for production usage.
`nkithr`, i.e., `nkith` revised, is recommended in its stead.

nkmax
-----

### Usage

mpirun -np <node count> `nkmax -- -n` <population count> [`-f` <dump file> `-g`]

`-n` <population size> : The count of numbers to be generated across *k* nodes.
`-f` <dump file> : The name of a binary file storing the generated population. Requires `-g`.
                   Default value is `./popl.out`.
`-g` : Debug mode. Distributively generates a dump file of the entire population using MP/IO.

### Notes

`nkmax` is an incredibly simple program. Each node generates its own local population.
As the population is generated, a marker is kept and updated for the maximum value seen heretofore.
When the generation is complete, all nodes send their local maxima to the root node.
The root node determines and reports the maximum of this set.

Time complexity: O(1)
Message complexity: O(k-1)

nkith
-----

### Usage

mpirun -np <node count> `nkith -- -n` <population count> `-i` <selection index>
                        [`-f` <dump file> `-g (--123|--321) -d` <duplicate count> ]
`-n` <population size> : The count of numbers to be generated across *k* nodes.
`-i` <selection index> : The index of the *i*th element to find.
`-f` <dump file> : The name of a binary file storing the generated population. Requires `-g`.
                   Default value is `./popl.out`.
`-g` : Debug mode. Distributively generates a dump file of the entire population using MP/IO.
`--123` : Generate a series of numbers ascending in step with *k* instead of a random set.
`--321` : Generate a series of numbers descending in step with *k* instead of a random set.
`-d` <duplicate count> : Duplicate each generated number this many times.
                         Implicit upper bound of *n*/*k* for randomly generated populations.

### Notes

`nkith` was a first attempt at solving the problem of finding the *i*th element by pre-sorting
the local population and doing cursor accounting instead of set manipulation. While it proved
reasonably quick for small values of *n*, the time cost ballooned as *n* scaled.

The program uses the C stdlib function `qsort` to achieve the local ordering.
Local pivot candidates are chosen deterministically by selecting the median of
the remaining live population. The root node keeps a round-robin counter and selects
the candidate from the node whose identifier corresponds to the counter.
The pivot-elect is broadcast to all nodes.

Starting from the *a priori* lower bound, each node counts until it reaches an element
greater than or equal to the pivot. It then counts the number of elements equal to
the pivot, if any exist in the local population. Two pieces of information are sent
back to the root node: the number of elements less than the pivot and the number less
than or equal to the pivot. The root node sums the respective quantities and broadcasts
the aggregated information back to all the nodes.

Given the cardinalities, each node updates its lower bound and upper bound cursors on
the live population. The process continues until the global count of elements lower
than the current pivot less than *i* AND the count of elements less than or equal to
the current pivot is at least *i*.

Provided that the system implementation of qsort is an in-place sort, this program
is conservative on memory since there would be only ever one copy of the local
population modulo the cursor variables required to perform the sort. However,
the in-place precondition is not guaranteed:
http://www.gnu.org/software/libc/manual/html_mono/libc.html#Array-Sort-Function

Note that the time bound is based on the nomenclature of the function and its
original implementation. It may in fact be worse depending on the system.

Time complexity: O([n/k] log [n/k]) [presort] + O(log n) [rounds]
Message complexity: O(k log n)

nkithr
------

### Usage

mpirun -np <node count> `nkithr -- -n` <population count> `-i` <selection index>
                        [`-f` <dump file> `-g (--123|--321) -d` <duplicate count> ]
`-n` <population size> : The count of numbers to be generated across *k* nodes.
`-i` <selection index> : The index of the *i*th element to find.
`-f` <dump file> : The name of a binary file storing the generated population. Requires `-g`.
                   Default value is `./popl.out`.
`-g` : Debug mode. Distributively generates a dump file of the entire population using MP/IO.
`--123` : Generate a series of numbers ascending in step with *k* instead of a random set.
`--321` : Generate a series of numbers descending in step with *k* instead of a random set.
`-d` <duplicate count> : Duplicate each generated number this many times.
                         Implicit upper bound of *n*/*k* for randomly generated populations.
`-D` : Deterministic pivot selection. Pivot selection is random by default.

`nkithr`, or `nkith` *r*evised (originally `nkith` *r*andom), is a proper
implementation of the distributed Random-Select algorithm described in
https://sites.google.com/site/gopalpandurangan/algbook.pdf?attredirects=0.

Instead of pre-sorting the locally generated populations, as `nkith` does,
the data movement in `nkithr` is limited to binning the live local population
into *S1*, the set of numbers below the pivot, or *S2*, the set of
numbers above the pivot. The set of numbers equal to the pivot is simply counted.

By default, the pivot candidates are selected at each node according to a uniform
distribution over the indices of the live local population. With the `-D`
command-line flag, the choice becomes deterministic, namely the index midpoint of
the live local population.

The candidates as well as the live local population size are sent to the root node
for final selection. The root node determines the proportional fractions residing
at each node with respect to the total live population and creates a discrete
distribution over [0,*k*) with the previously mentioned fractional weights.
In the random case, the pivot is chosen from the candidates according to this
distribution. In the deterministic case, the round-robin method from `nkith` is
reused with a slight modification: nodes are skipped if their corresponding
weight is zero (0), and the counter is incremented until a non-zero weight node
is found. This prevents wasting iterations on nodes which have no truly viable
pivot candidates to supply.

When each node receives the pivot-elect, it walks from the lower bound (inclusive)
to the upper bound (exclusive) on the population array *M*. The bounds start
on the interval [0,*m*) where *m* is the local population. Values less than the
pivot are placed into the set *S1*, while values greater than the pivot are
placed into *S2*. Values equal to the pivot are counted. The storage of
*S1* and *S2* is implemented by an array *Q* of size *m* where elements
belonging to *S1* are inserted with an incrementing cursor at the start
of the array while elements belonging to *S2* are inserted with a decrementing
cursor at the end of the array. When the upper bound is reached, the array
is filled with the pivot for |*P*| cells after the last element of *S1*
where *P* is the set of elements equal to the pivot. Thus, *Q* becomes a
complete semi-sorted copy of *M* after the first iteration.

The local cardinalities |*S1*| and |*S1*| + |*P*| are sent to the root node from every
node, and the total cardinalities are aggregated and rebroadcast to all nodes.
The live subspace from the upper bound to the lower bound is copied from *Q* to *M*
to preserve the semi-sort for the next round.

If the global cardinality |*S1*| > *i*, then the upper bound is set at |*S1*|.
If |*S1*| + |*P*| < i, then the lower bound is set at |*S1*| + |*P*|.

While this requires double the natural memory required to store a given population,
it saves tremendously on data movement, and as the live subspace shrinks, the algorithm
becomes increasingly cache-friendly. It was for the sake of avoiding memory and cache churn
that the decision was made not to free the old *Q* and reallocate a smaller *Q* in every
iteration. Even if that path had been chosen, the initial doubling spike in the first
iteration would have occurred regardless. Time was favored at the expense of local memory.

The termination condition is the same as in `nkith`: the global count of elements lower
than the current pivot must be less than *i* AND the count of elements less than or equal to
the current pivot must be at least *i*.

It should be noted that the time required for population generation is a substantial portion
of the execution (wall-clock) time. It is suspected that this is also the case for `nkith`,
but since the generation methodology is the same, the dominating term in the time equation
there is most likely the local pre-sort. However, `nkith` only provides measurements for
total execution time.

It is not entirely accurate to call the deterministic pivot selection option fully deterministic.
Since the numbers themselves are generated randomly (unless specified otherwise at the command
line) and placed in no certain order in the array *M*, choosing the index midpoint has a certain
element of randomness to it.

Time complexity: O(log n) [rounds]
Message complexity: O(k log n)

popldump
--------

### Usage

[mpirun -np <node count>] `popldump -- -n` <population count> `-i` <selection index>
                        [`-f` <dump file> ]
`-n` <population size> : The count of numbers to print.
`-f` <dump file> : The name of a binary file storing the generated population. Requires `-g`.
                   Default value is `./popl.out`.

### Notes

`popldump` is simply the inversion of the `-f` and `-g` flags for the `nk`* family of programs.
It prints out each member of the population stored in a given dump file on a separate line.
This is useful for checking the correctness of the output of `nk`*. Consider the following example:

`mpirun -np 2 nkithr -- -n 16 -i 8 -g -f test.out; ./popldump -n 16 -f test.out | sort -n | head -8 | tail -1`

The number under the `p` column should match the single number printed underneath, e.g.,

`n,k,i,p,m,r,t,T`
`16,2,8,1781118644,25,6,172421,5404230`
`1781118644`

Note that invocation with `mpirun` is not necessary. This is because all the work is done
at the root node, namely node 0. There is no guarantee in the distributive case that stdout
will be line-atomic, which is why the MP/IO debug dumps were built into the `nk`* family of
programs in the first place.

Users are recommended to be mindful of the time cost of doing dumps of large populations.
The execution of the writer may not take long, but a reader like

`./popldump -n $[2 ** 28] | sort -n | head -$[2 ** 27] | tail -1`

will provide an excellent opportunity to get up and make a cup of coffee, walk the dog,
do taxes, finish a graduate degree, etc.

Known Issues
------------

* `popldump` does not guard against exceeding a file's size with `-n`.
* `nkith`* do not guard against `-i` > `-n`.
* `nkith`* and `nkmax` do not have a means to read a dump file and distribute to a set of nodes.

Contact
-------

Please send comments, questions, and bug reports to Michael Yantosca via e-mail at
mike@archivarius.net.
