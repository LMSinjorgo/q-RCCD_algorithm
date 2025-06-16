# $q$-RCCD algorithm
## By A. Ghaffari-Hadigheh, L. Sinjorgo, and R. Sotirov (2024)

Implementation of the $q$-RCCD algorithm, for solving problems of the form

$$\min  f(x)$$ 
$$\text{subject to: } a'x = b, \text{ } l \leq x \leq u,$$ 

where $f$ is a differentiable, possibly non-convex, function. In each iteration, we randomly choose $q$ coordinates of $x$ (uniform probability).
Then, based on a convex approximation of $f$, we solve a simple quadratic programme on the selected variables, and update $x$ accordingly.

This repository contains implementations of the $q$-RCCD algorithm, applied to the densest $k$-subgraph ($\text{D}k\text{S}$) problem, and the eigenvalue complementarity (EiC) problem.

More details and a convergence analysis are provided in our paper, available [here.](https://link.springer.com/article/10.1007/s10898-024-01429-6)
```
@article{ghaffari2024convergence,
  title={On convergence of aq-random coordinate constrained algorithm for non-convex problems},
  author={Ghaffari-Hadigheh, Alireza and Sinjorgo, Lennart and Sotirov, Renata},
  journal={Journal of Global Optimization},
  volume={90},
  number={4},
  pages={843--868},
  year={2024},
  publisher={Springer}
}
```
