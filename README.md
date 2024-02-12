# SVD Routine

### Date: December 9, 2017

### Author: Hong Hun

## Abstract
Singular Vector Decomposition (SVD) facilitates the diagonalization of a matrix A, comprising the Left Singular Vector (U), Right Singular Vector (V), and Singular Values (S). Essentially, A is reconstructed through the matrix multiplication of U, S, and the transpose of V. Singular vectors in U and V correspond to the eigenvectors of AA^T and A^TA, respectively, with S being a diagonal matrix containing the square roots of the eigenvalues from AA^T and A^TA (Banerjee et al., 2014). The Jacobi method, utilizing rotation matrices for diagonalization, is adopted for eigenvector computation (Drmac, 1999).

## Background
SVD computation is supported across various software packages that offer optimized algorithms for enhanced efficiency. Despite their power and speed, these packages often come with complex codebases that pose challenges in understanding and modification. To address these issues, this work presents a straightforward implementation in Fortran, devoid of optimization complexities, to aid in grasping statistical mathematics and coding.

## Method
The Jacobi method, an iterative algorithm, is utilized for eigenvector computation, transforming the m x n matrix A via rotational matrices P. This process gradually reduces off-diagonal elements, ultimately diagonalizing A to represent its eigenvalues.

## Code Implementation
Implemented in Fortran 90, the coding strategy employs sine, cosine, and absolute value functions, eschewing the need for additional built-in functions. The `svdc` main program orchestrates the input of matrix A and oversees the execution of SVD computation subroutines, including eigenvector calculations and matrix manipulations.

## Results
The Jacobi method's computational efficiency is benchmarked against an alternative SVD routine, highlighting its increased computational demand with larger matrix dimensions. Despite this, the method's simplicity and educational value in understanding algorithms are underscored.

## References
- Banerjee, Sudipto; Roy, Anindya. Linear Algebra and Matrix Analysis for Statistics. Chapman and Hall/CRC, 1st edition, 2014.
- Drmac, Z. 'A posteriori computation of the singular vectors in a preconditioned Jacobi SVD algorithm.' IMA Journal of Numerical Analysis, vol. 19, no. 2, April 1999, pp. 191–213.
