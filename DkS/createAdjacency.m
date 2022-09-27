function A = createAdjacency(n, e)
% creates a random adjacency matrix A, corresponding to a graph on n
% vertices, with e edges

matrixIdx = find(triu(ones(n,n),1));
matrixIdx = randsample(matrixIdx,e,false);
density = 2*e/n^2;


if density < 0.25
    % create sparse
    A = sparse(n,n);
    A(matrixIdx) = 1;
    A = A+A';
else
    A = zeros(n,n);
    A(matrixIdx) = 1;
    A = A+A';
end

end

