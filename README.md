# $q$-RCCD algorithm
## By A. Ghaffari-Hadigheh, L. Sinjorgo, and R. Sotirov (2023)

Implementation of the $q$-RCCD algorithm, for solving problems of the form
$$\min  f(x)$$ 

$$\text{subject to: } a'x = b, \text{ } l \leq x \leq u,$$ 
where $f$ is a differentiable, possibly non-convex, function. In each iteration, we randomly choose $q$ coordinates of $x$ (uniform probability).
Then, based on a convex approximation of $f$, we solve a simple quadratic programme on the selected variables, and update $x$ accordingly.

This repository contains implementations of the $q$-RCCD algorithm, applied to the densest $k$-subgraph (\text{D}k\text{S}) problem, and the eigenvalue complementarity (EiC) problem.

More details and a convergence analysis are provided in our paper.


```
@article{ghaffari2023convergence,
  title={On convergence of a $q$-random coordinate constrained algorithm for non-convex problems},
  author={Ghaffari-Hadigheh, Alireza and Sinjorgo, Lennart and Sotirov, Renata},
  journal={arXiv preprint arXiv:2210.09665},
  year={2023}
}
```
