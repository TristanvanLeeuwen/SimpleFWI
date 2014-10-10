SimpleFWI
=========

This code provides the basic building blocks to test optimization algorithms on seismic inverse problems.
The canonical seismic waveform inversion problem is given by

$$\min_{m} \sum_{i} \frac{1}{2}||P^TA(m)^{-1}q_i - d_i||_2^2 + \frac{\alpha}{2}||Lm||_2^2,$$

where $A(m)$ is the discretized Helmholtz operator $\omega^2 m + \nabla^2$ with absorbing boundary conditions, 
$P$ is the sampling operator, $q_i$ are the source vectors, $d_i$ are the observed data and $L$ is the discretized 
$\nabla$.

The main function is *misfit.m*, which takes as input the medium parameters $m$ and the data $d_i$ (and definitions of $P$, $q_i$, $\omega$ etc.) and returns the misfit value, its gradient and a function handle to compute the action of the Gauss-Newton hessian.

For an example, see */examples/marm.m*, which uses a simple BB iteration to solve the optimization problem. 
Replace this with your favourite optimization algorithm and you're ready to roll!

Feel free to contact me with any questions, suggestions, etc.

Tristan van Leeuwen - T.vanLeeuwen@uu.nl
 
