SimpleFWI
=========

This code provides the basic building blocks to test optimization algorithms on seismic inverse problems.
The canonical seismic waveform inversion problem is given by

$$\min_{m} \sum_{i} \frac{1}{2}||P^TA(m)^{-1}q_i - d_i||_2^2 + \frac{\alpha}{2}||Lm||_2^2,$$

where $A(m)$ is the discretized Helmholtz operator $\omega^2 m + \nabla^2$ with absorbing boundary conditions, 
$P$ is the sampling operator, $q_i$ are the source vectors, $d_i$ are the observed data and $L$ is the discretized 
$\nabla$.

The main function is |misfit.m|
