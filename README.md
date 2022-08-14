# Triangular-SUB

Triangular systems of equations are of the form

$$ Ax = b $$

where $A \in \mathbb{R}^{n \times n}$ a known *upper* or *lower triangular* matrix, $b \in \mathbb{R}^{n}$ a known vector, and $x \in \mathbb{R}^{n}$ the vector we want to solve the system for.

## Motivation

Most introductory courses on numerical analysis introduce two algorithms for solving triangular systems of equations. Usually these algorithms are implemented using two nested (uncompiled) for-loops. A significant performance improvement can be made by rewriting the inner of these two loops in terms of a vector operation (i.e., a compiled for-loop).

Commonly, the standard forward and backward substitution algorithms (for solving lower and upper triangular systems respectively) are treated as two completely different entities. Below, my aim is to demonstrate that a vectorized implementation of these algorithms happens to yield the exact same code up to one little difference in the indexing.

## Mathematical background to the implementation

Let me write a lower triangular system in an explicit (pythonic) notation:

$$ \begin{bmatrix} A[0,0] & & & & & \\\ A[1,0] & A[1,1] & & & & \\\ \vdots & \vdots & \ddots & & & \\\ A[i,0] & A[i,1] & \dots &  A[i,i] & & \\\ \vdots & \vdots & & \vdots & & \ddots \end{bmatrix} \begin{bmatrix} x[0] \\\ x[1] \\\ \vdots \\\ x[i] \\\ \vdots \end{bmatrix} = \begin{bmatrix} b[0] \\\ b[1] \\\ \vdots \\\ b[i] \\\ \vdots \end{bmatrix} $$

Suppose $x[0], x[1], \dots, x[i-1]$ are already known. By explicit multiplication of the $i$-th row in $A$ with $x$, we see that

$$ b[i] = \sum_{j=0}^{i} A[i,j] x[j] $$

Introduce a helper vector

$$ \chi^{(i)} = \begin{bmatrix} x[0] \\\ x[1] \\\ \vdots \\\ x[i-1] \\\ 0 \\\ \vdots \\\ 0 \end{bmatrix} $$

which allows us to rewrite the above equation in terms of a dot product (denoted by $\langle, \rangle$):

$$ b[i] = \sum_{j=0}^{i} A[i,j] x[j] = \sum_{j=0}^{i-1} A[i,j] x[j] + A[i,i] x[i] = \langle A[i,:], \chi^{(i)} \rangle + A[i,i] x[i] $$ 

Here, I have employed the *numpy* indexing notation $A[i,:]$ to denote the $i$-th row of $A$.

Some rearranging then yields

$$ x[i] = (b[i] - \langle A[i,:], \chi^{(i)} \rangle) / A[i,i] $$

Hence, starting with $\chi^{(0)}$ (which by definition is just the zero vector), we can iteratively "fill" in this helper vector's components to finally obtain the solution $x$ to the system. But instead of explicitly working with a helper vector, we can just start with the vector $x$ directly by first setting all its components to zero and then iteratively filling in one component after the other according to the above equation.

An almost identical study can be conducted for the case of an upper triangular matrix, where the only difference is that one starts from the last component of $x$ and then iteratively fills in the previous components in the inverted order.

## Algorithms

### Forward substitution

Inputs: Lower triangular matrix $A$, vector $b$
1. Let $x$ be the zero vector
2. For $i = 0, 1, \dots, n-1$ let $x[i] = (b[i] - \langle A[i,:], x \rangle) / A[i,i]$

### Backward substitution

Inputs: Lower triangular matrix $A$, vector $b$
1. Let $x$ be the zero vector
2. For $i = n-1, n-2, \dots, 0$ let $x[i] = (b[i] - \langle A[i,:], x \rangle) / A[i,i]$

One sees immediately that the two algorithms are basically the same, except that the order of the loop is inverted. This motivated me to combine the algorithms.
