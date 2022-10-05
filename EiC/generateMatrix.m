function A = generateMatrix(n,d)
% Generates a matrix suitable for EiC.
% output A is a symmetric matrix of size n, with density d and positive
% diagonal

% recompute density because we add positive diagonal
% elements after.
dPrime = (d-1/n);
A = sprand(n,n,dPrime/2);
A = A'+A;
A = A+spdiags(0.01+abs(normrnd(0,1,n,1)),0,n,n);
end