function A = generateMatrix(n,d)
% Generates a matrix suitable for EiCP
% about p nonzeros per row

% recompute density because we add positive diagonal
% elements after.

A = abs(sprandsym(n,d));

% ensure positive diagonal
A = A+spdiags(0.01+abs(normrnd(0,1,n,1)),0,n,n);

end

