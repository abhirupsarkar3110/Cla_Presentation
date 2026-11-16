# Algorithm Results: Inverse Norm of Tridiagonal Matrices

## Overview

This document tracks the results and performance of the matrix norm algorithms `nrminv_final1` and `nrminv_final2`. These algorithms are designed to handle zero/infinite pivots and block pivoting for robust computation of the inverse norm.

## Test Case Configuration

For our initial testing, we used a $5 \times 5$ symmetric positive-definite tridiagonal matrix $A$ (a 1D discrete Laplacian), structured as follows:

* **Main diagonal (**$a$**):** $2.0$
* **Sub-diagonal (**$b$**):** $-1.0$
* **Super-diagonal (**$c$**):** $-1.0$

Mathematically, the matrix $A$ is represented as:

$$
A = \begin{pmatrix} 2 & -1 & 0 & 0 & 0 \\ -1 & 2 & -1 & 0 & 0 \\ 0 & -1 & 2 & -1 & 0 \\ 0 & 0 & -1 & 2 & -1 \\ 0 & 0 & 0 & -1 & 2 \end{pmatrix}
$$

## Results

When executed against the test matrix, the calculated infinity-norm of the inverse is as follows:

### 1. `nrminv_final1.py` (NaN-free, handles zero/infinite pivots)
* **Status:** Success
* **Calculated Norm (**$\gamma$**):** `4.5`

### 2. `nrminv_final2.py` (2x2 block pivoting, no division by tiny pivots)
* **Status:** Success
* **Calculated Norm (**$\gamma$**):** `4.5`

## Conclusion
Both algorithms successfully handled the inversion and correctly calculated the infinity-norm of the matrix's inverse. The result of $4.5$ matches the exact mathematical expectation for this specific matrix configuration.