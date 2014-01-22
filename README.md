Topics in Randomized Numerical Linear Algebra
=============================================

My PhD thesis, defended in May 2013 at the California Institute of Technology under the advisement of Professor Joel Tropp. 

Abstract
--------

This thesis studies three classes of randomized numerical linear algebra algorithms, namely: 
* randomized matrix sparsification algorithms, 
* low-rank approximation algorithms that use randomized unitary transformations, and 
* low-rank approximation algorithms for positive-semidefinite (PSD) matrices.

Randomized matrix sparsification algorithms set randomly chosen entries of the input matrix to zero. When the approximant is substituted for the original matrix in computations, its sparsity allows one to employ faster sparsity-exploiting algorithms. This thesis contributes bounds on the approximation error of nonuniform randomized sparsification schemes, measured in the spectral norm and two NP-hard norms that are of interest in computational graph theory and subset selection applications.

Low-rank approximations based on randomized unitary transformations have several desirable properties: they have low communication costs, are amenable to parallel implementation, and exploit the existence of fast transform algorithms. This thesis investigates the tradeoff between the accuracy and cost of generating such approximations. State-of-the-art spectral and Frobenius-norm error bounds are provided.

The last class of algorithms considered are SPSD “sketching” algorithms. Such sketches can be computed faster than approximations based on projecting onto mixtures of the columns of the matrix. The performance of several such sketching schemes is empirically evaluated using a suite of canonical matrices drawn from machine learning and data analysis applications, and a framework is developed for establishing theoretical error bounds.

In addition to studying these algorithms, this thesis extends the Matrix Laplace Transform framework to derive Chernoff and Bernstein inequalities that apply to all the eigenvalues of certain classes of random matrices. These inequalities are used to investigate the behavior of the singular values of a matrix under random sampling, and to derive convergence rates for each individual eigenvalue of a sample covariance matrix.
