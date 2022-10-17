# qRCCA-algorithm
By A. Ghaffari-Hadigheh, L. Sinjorgo, and R. Sotirov
Implementation of the $q$ RCCA algorithm, for solving problems of the form
$$\max  f(x)$$ 

$$\text{subject to: } a'x = b, \text{ } l \leq x \leq u,$$ 
where $f$ is differentiable, possibly non-convex, function. In each iteration, we randomly choose $q$ coordinates of $x$ (uniform probability).
Then, based on a convex approximation of $f$, we solve a easy quadratic programme on the selected variables, and update $x$ accordingly.

For a convergence analysis and more details, see our paper.
