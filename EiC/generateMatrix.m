function A = generateMatrix(n,d)
% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
%
% Generates a matrix suitable for EiC.
% output A is a symmetric matrix of size n, with density d and positive
% diagonal

% recompute density because we add positive diagonal
% elements after.
dPrime = (d-1/n);
A = sprand(n,n,dPrime/2);
A = A'+A;
A = A+spdiags(0.001+abs(normrnd(0,1,n,1)),0,n,n);
end