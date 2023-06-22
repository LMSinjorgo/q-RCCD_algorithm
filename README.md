# q-RCCD algorithm
## By A. Ghaffari-Hadigheh, L. Sinjorgo, and R. Sotirov (2023)

Implementation of the $q$-RCCD algorithm, for solving problems of the form
$$\min  f(x)$$ 

$$\text{subject to: } a'x = b, \text{ } l \leq x \leq u,$$ 
where $f$ is differentiable, possibly non-convex, function. In each iteration, we randomly choose $q$ coordinates of $x$ (uniform probability).
Then, based on a convex approximation of $f$, we solve a easy quadratic programme on the selected variables, and update $x$ accordingly.

More details and a convergence analysis are provided in our paper.


% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
