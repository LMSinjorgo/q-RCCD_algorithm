function [objValue, x,iterCount] = qRCCD_DkS(A,q,k,maxTime,maxIter,x)
% official Code for the paper:
% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
%
% Implementation of the q-RCCD algorithm for the DkS problem
%
% INPUTS:
% A: symmetric adjacency matrix of graph
% q: number of coordinates to update simultaneously
% k: number of vertices to select (appears in the constraint sum(x) = k )
% maxTime: maximum running time (seconds)
% maxIter: maximum number of iterations
% x: feasible starting vector (Not necessary, can be left empty).
%
% Outputs are obvious.
%
% DkS:
% maximize: x'Ax,
% subject to: sum(x) = k, 0 \leq x \leq 1. Here, A is the adjaency matrix
% of a graph
%
% In each step, we solve a projection Problem:
% u = proj(nabla f / L + x),
% using the code fullProj (see references therein)

% input checking
givenX = false;
if nargin == 6
    if ~isempty(x)
        givenX = true;
    end
end

if givenX
    if sum(x) ~= k || any(x > 1) || any(x < 0)
        error("The provided x vector is not feasible!");
    end
    if size(x,1) == 1
        x = x';
    end
end

% initial feasible starting point (if no starting point is given by the
% user)
n = size(A,2);
if ~givenX
    x = (k/n) * ones(n,1);
end

if ~issymmetric(A)
    error("matrix A is not symmetric");
end

if ~islogical(A)
    error("Given adjacency matrix A is not logical! Please make this a logical matrix for speed!")
end

if numel(x) ~= n
    error("Vector x is of the wrong length!")
end

if q > n
    error("Value of q > n! q should be n or smaller.")
end

if isempty(maxIter)
    maxIter = Inf;
end
if isempty(maxTime)
    maxTime = Inf;
end

if maxIter < 0
    error("maxIter is negative!")
end

if maxTime < 0
    error("maxTime is negative!")
end

if isinf(maxIter) && isinf(maxTime)
    error("No valid stopping criteria given!")
end


% start timer
tic

% store value of A*x for speed
linA = A*x;

iterCount= 0;

zeroVector = zeros(q,1);

while toc < maxTime && iterCount < maxIter
    J = uint32(randperm(n,q));
    xJ = x(J);

    % compute parameters for the projection problem
    sum_xJ = sum(xJ);

    % if sum_xJ = 0, then the projection must be the all 0's vector
    % so we can skip this iteration
    if sum_xJ == 0
        iterCount = iterCount+1;
        continue;
    else
        % compute L as || A_J ||_1
        colA = A(:,J);
        L = 2*full(max(sum(colA(J,:),2)));
        L = max(L,10^(-5));
        % compute projection problem
        % Projection algorithm used from:
            % Wang, Weiran, and Canyi Lu. 
            % "Projection onto the capped simplex." 
            % arXiv preprint arXiv:1503.01002 (2015).
            % (we have slightly adapted the code for correct handling of
            % edge cases, and remove unnecessary operations for our end
            % goals)
            % code taken from: 
            % https://github.com/canyilu/Projection-onto-the-capped-simplex
        uJ=fullProj((2/L)*linA(J) + xJ,sum_xJ,q,zeroVector);
    end

    % compute dJ and the next iterate of x
    x(J) = uJ;
    dJ = uJ - xJ;

    % update the variable tracking Ax
    linA = linA + colA(:,dJ ~= 0)*dJ(dJ ~= 0);

    iterCount = iterCount+1;
end

% Output checking
if abs(sum(x)-k) >= 0.000001
    error("Vector x not feasible!")
end

objValue = x'*A*x;
end