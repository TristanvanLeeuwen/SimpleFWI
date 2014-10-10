SimpleFWI
=========

This code provides the basic building blocks to test optimization algorithms on seismic inverse problems.
The canonical seismic waveform inversion problem is given by

$$\min_{m} \sum_{i} \frac{1}{2}||P^TA(m)^{-1}q_i - d_i||_2^2$$

