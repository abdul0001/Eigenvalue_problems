# Eigenvalue_problems
We will design MATLAB code for solving Eigenvalue Problems


In vectorized_P1FEM_EVP we have implemented the Laplace eigenvalue problem using P1FEM and alsp implemented the ROM.

There are two main files main_unform.m where we have used uniformly distributed time steps for generating solution and in the 
snapshot matrix we have take all the solutions of 4,8,... iteratins. 

In (main_non_uniform.m) the non-uniform case we have generated Latin hypercube sampling to generate values between [0,10] 
and used them for solving the eigenvalue problem and used all the solutions as snapshots.
