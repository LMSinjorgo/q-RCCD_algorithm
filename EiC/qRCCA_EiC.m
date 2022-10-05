function [runningTime, objValue, x] = qRCCA_EiC(A,B,q,M)
% Implementation of the qRCCA algorithm for EiC.
% Lennart Sinjorgo & Renata Sotirov

% inputs:
% A,B: symmetric sparse matrices with positive diagonal
% Density of these matrices should be low, to ensure a fast run time

% q: number of coordinates to update simultaneously

% M: Number of iterations to run (no other stopping criteria)


% EiC:
% minimize: log(x'Ax) - log(x'Bx),
% subject to: sum(x) = 1, x \geq 0

% In each step, we solve a QP for the update:
% max f(u)  = NablaF (u - xJ) - L/2 * || u - xJ ||^2
%           = NablaF *u - L/2 * (u'u - 2u'x)      + Constants
%           = - u' *(L/2)* u + (NablaF + L*x)'*u  + Constants
%           = 0.5 u'*( -L*I )*u + c'u                  + Constants
% s.t. aJ * uJ = aJ * xJ

tic
n = size(A,2);

% starting Value
x = (1/n) * ones(1,n);

% initialize parameters
normA = full(sum(A,1));
normB = full(sum(B,1));

% quadFormA = x'*A*x;
% quadFormB = x'*B*x;
quadFormA = sum(normA) / (n^2); 
quadFormB = sum(normB) / (n^2); 

onesColVector = (1:q)';

m = 1;

% Matrix G contains the first q rows of matrix:
% [             1 ] ^(-1)
% [     I       1 ]
% [             1 ]
% [ 1   1    1  0 ]
G = [eye(q) - (1/q)*ones(q), (1/q)*ones(q,1)];

% compute 1-norms of matrix A and B
normA = max(normA);
normB = max(normB);
M = M+1;

% store the diagonals for faster indexing
diagA = full(diag(A));
diagB = full(diag(B));

% linIdxC and linIdxR are used for linear indexing of matrices A and B
linIdxC = repelem(1:(q-1),(q-1):-1:1);
linIdxR = [];
for p = 2:q 
    linIdxR = [linIdxR, p:q];
end

% sIdx is used for linear indexing of dJ*dJ' (the cross terms)
sIdx = sub2ind([q q],linIdxR,linIdxC)';

while m <  M
    % update randomly q coordinates
    J = randperm(n,q);
    xJ = x(J);
    
    % use function ind2sub
    % without all the input checking by MatLAB
    idx = J(linIdxR) + (J(linIdxC)-1).*n;
    
    % compute Lipschitz constant
    % using || A_J || <= || A ||
    % and   || B_J || <= || B ||
    L = 2*( normA /quadFormA + normB / quadFormB );

    % store intermediate values
    z_A = x*A(:,J);
    z_B = x*B(:,J);

    % compute partial derivatives
    NablaF = 2 * ( z_A / quadFormA - z_B / quadFormB );

    % compute parameters
    sum_xJ = sum(xJ);
    linearCoeff = NablaF/L + xJ;

    % project onto part of the feasible set
    uJ = G * [linearCoeff, sum_xJ]';
    
    % if not feasible, project again.
    % Projection code by
    % Condat, Laurent. 
    % "Fast projection onto the simplex and the l1 Ball." 
    % Mathematical Programming 158.1 (2016): 575-585
    if any(uJ < 0)
        uJ = max(uJ-max((cumsum(sort(uJ,1,'descend'),1)-sum_xJ)./onesColVector),0);
    end

    % compute dJ and the next iterate of x
    uJ = uJ';
    dJ = uJ - xJ;
    x(J) = uJ;

    % update the variables tracking x'*A*x and x'*B*x
    dSquared = dJ.^2;

    % find the off diagonal elements in A_J
    offDiagA = A(idx);
    % if there are none, then we only need the diagonal entries of A
    if nnz(offDiagA) == 0
        quadFormA = quadFormA + 2*z_A*dJ'+dSquared * diagA(J);
        outerProdUpdate = 0;
    else
        outerProd = dJ' * dJ;
        quadFormA = quadFormA + 2*z_A*dJ'+dSquared * diagA(J) + 2*offDiagA*outerProd(sIdx);    
        outerProdUpdate = 1;
    end

    offDiagB = B(idx);
    if nnz(offDiagB) == 0
        quadFormB = quadFormB + 2*z_B*dJ'+dSquared * diagB(J);
    else
        if outerProdUpdate == 0
            outerProd = dJ'*dJ;
        end
        quadFormB = quadFormB + 2*z_B*dJ'+dSquared * diagB(J) + 2*offDiagB*outerProd(sIdx);    
    end
    
    
    m = m+1;
end
runningTime = toc;

% output objective value
objValue = log(quadFormA/quadFormB);
end