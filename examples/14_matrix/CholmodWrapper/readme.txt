We offer some examples to solve a linear least square problem Ax=b by CholmodWrapper (http://www.mpi-inf.mpg.de/~hongbofu/software.html#cholmod_wrapper), here.

CholmodWrapper is under a LGPL license.

There is a thorough performance comparison between CHOLMOD and other sparse solvers (e.g. TAUCS) by [Gould et al. 2005]
(ftp://ftp.numerical.rl.ac.uk/pub/reports/ghsNAGIR20051r1.pdf).
CHOLMOD and TAUCS are among the fastest direct solvers. Specifically,CHOLMOD is the fastest for 42 of the 87 matrices, which they used for
testing. Its running time is either fastest or within 10% of the fastest for 73 out of 87 matrices.

Sorry, I didn't. But actually there is one SIGGRAPH paper 2006 in which the authors did a performance comparison between TAUCS and
CHOLMOD for the mesh deformation application. If you're interested, you can download the paper at http://www-sal.cs.uiuc.edu/~yyz/publication/multigrid_mesh.pdf. From the data they provided, I would say TAUCS and CHOLMOD is comparable,
at least.

One unique feature of CHOLMOD I like is its support of Cholesky factorization modification (update/downdate), which is very useful for
problems with changing constraints, e.g. interactive mesh deformation.


The first exampe is:
	        [ 1	1 ]
	 A =    [ 1	-1]
	        [ 1	2 ]
   
         b=    [ 2 ]
	       [ 0 ]
               [ 4 ]

Using Matlab (least_square_examples.m), we get the solution:

                [ 1.5]
	 x =    [ 1.5]

Using the c++ (example)
