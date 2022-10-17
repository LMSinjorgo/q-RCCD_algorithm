function [runningTime, objValue, x] = qRCCA_DkS(A,q,k,M,x,eps,maxStag)
% Implementation of the qRCCA algorithm for
% the continuous relaxation of DkS (densest k subgraph problem).
% Lennart Sinjorgo & Renata Sotirov, 2022
%
% INPUTS:
% A: symmetric adjacency matrix of graph
% q: number of coordinates to update simultaneously
% k: number of vertices to select (appears in the constraint sum(x) = k )
% M: maximum number of iterations (no other stopping criteria)
% x: feasible starting vector (Not necessary, can be left empty).
% eps and maxStag: if sum(abs(x_k - x_{k+1}) < eps,
%                  counter = counter + 1. If counter > maxStag, stop iterating
%
% Outputs are obvious.
%
% DkS:
% minimize: x'Ax,
% subject to: sum(x) = k, 0 \leq x \leq 1. Here, A is the adjaency matrix
% of a graph
%
%
%
% In each step, we solve a QP for the update:
% max f(u)  = NablaF (u - xJ) - L/2 * || u - xJ ||^2
%           = NablaF *u - L/2 * (u'u - 2u'x)      + Constants
%           = - u' *(L/2)* u + (NablaF + L*x)'*u  + Constants
%           = 0.5 u'*( -L*I )*u + c'u                  + Constants
% s.t. aJ * uJ = aJ * xJ


% input checking
if ~isempty(x)
    if sum(x) ~= 1 || any(x > 1) || any(x < 0)
        error("The provided x vector is not feasible!");
    end
end

if ~issymmetric(A)
    error("matrix A is not symmetric");
end

tic
n = size(A,2);

% initialize feasible starting point (if no starting point is given by the
% user)
if isempty(x)
    x = (k/n) * ones(n,1);
end

% store value of Ax for fast indexing
linA = A*x;

m = 0;

% compute Lipschitz constant
% i.e., || A ||_1
L = full(max(sum(A,2)));

% Matrix G contains the first q rows of matrix:
% [             1 ] ^(-1)
% [     I       1 ]
% [             1 ]
% [ 1   1    1  0 ]
G = [eye(q) - (1/q)*ones(q), (1/q)*ones(q,1)];

% store zeroVector for speed
zeroVector = zeros(q,1);
stagnationCounter = 1;
while m < M
    J = randperm(n,q);
    xJ = x(J);

    % compute partial derivatives
    NablaF = 2 * linA(J,:);

    % compute parameters for the Quadratic programme
    sum_xJ = sum(xJ);

    % if sum_xJ = 0, then the projection must be the all 0's vector
    if sum_xJ == 0
        uJ = zeros(q,1);
    else
        linearCoeff = NablaF/L + xJ;

        % compute first projection onto sum(u) = sum(xJ)
        uJ = G * [linearCoeff; sum_xJ];
        
        % if this projection is infeasible, do second projection onto set
        % sum(u) = sum(xJ), 0 <= u <= 1
        if any(uJ < 0) || any(uJ > 1)
            % Projection algorithm used from
            % Wang, Weiran, and Canyi Lu. 
            % "Projection onto the capped simplex." 
            % arXiv preprint arXiv:1503.01002 (2015).
            % (we have slightly adapted the code for correct handling of
            % edge cases, and remove unnecessary operations for our end
            % goals)
            % code taken from: 
            % https://github.com/canyilu/Projection-onto-the-capped-simplex
            % We have also implemented the MEX version of function
            % fullProj, but this led to no difference in speed.
            uJ = fullProj(uJ,sum_xJ,q,zeroVector);
        end
    end

    % compute dJ and the next iterate of x
    dJ = uJ - xJ;
    
    x(J) = uJ;

    % update the variable tracking Ax
    linA = linA + A(:,J)*dJ;

    m = m+1;

    % check stopping criteria
    if sum(abs(dJ)) < eps
        stagnationCounter = stagnationCounter+1;
        if stagnationCounter > maxStag
            break;
        end
    end
end

runningTime = toc;

objValue = x'*A*x;
end
