function [objValue, x,iterCount] = qRCCD_EiC(A,B,q,maxTime,maxIter,x)
% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
%
% Implementation of the q-RCCD algorithm for the EiC problem
%
% Inputs:
% A,B: symmetric sparse matrices with positive diagonal
% q: number of coordinates to update simultaneously
% maxTime: maximum running time in seconds
% maxIter: maximum number of iterations
% x: feasible initial point, can be left empty, in which case the code
% will provide a starting point.
%
% EiC:
% max: log( (x'Ax) / (x'Bx) )
% subject to: sum(x) = 1, x \geq 0
%
% In each step, we solve a projection Problem:
% u = proj(nabla f / L + x),
% using the code provided by
% Condat, Laurent. 
% "Fast projection onto the simplex and the l1 Ball." 
% Mathematical Programming 158.1 (2016): 575-585    

% input checking
givenX = false;
if nargin == 6
    if ~isempty(x)
        givenX = true;
    end
end

if givenX
    if sum(x) ~= 1 || any(x > 1) || any(x < 0)
        error("The provided x vector is not feasible!");
    end
    if size(x,2) == 1
        x = x';
    end
end

% initial feasible starting point (if no starting point is given by the
% user)
n = size(A,2);
if ~givenX
    x = (1/n) * ones(1,n);
end

if ~issymmetric(A) || ~issymmetric(B)
    error("Input matrices not symmetric");
end

if any(diag(A) <= 0) || any(diag(B) <= 0)
    error("matrices contain nonpositive diagonal!")
end

if any(A(:) < 0) || any(B(:) < 0)
    error("matrices contain negative elements!")
end

if numel(x) ~= n
    error("vector x is of the wrong length")
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

% initialize parameters
quadFormA = x*A*x';
quadFormB = x*B*x';

onesRowVectorInv = (1:q).^(-1);
iterCount = 1;


while toc <  maxTime && iterCount < maxIter
    % update randomly q coordinates
    J = uint64(randperm(n,q));
    xJ = x(J);

    sumXJ = sum(xJ);
    if sumXJ == 0
        iterCount=iterCount+1;
        continue;
    end

    % store intermediate values
    % (index columns, not rows, for faster acces)
    A_J = A(:,J);
    z_A = x*A_J;

    B_J = B(:,J);
    z_B = x*B_J;

     % compute || A_J || and  || B_J ||
    subA = A_J(J,:); 
    LA = full(max(sum(subA,1)));
    subB = B_J(J,:); 
    LB = full(max(sum(subB,1)));

    % compute partial derivatives
    halfNablaF = ( z_A / quadFormA - z_B / quadFormB );

    % Compute Lipschitz constant
    halfL = ( LA / quadFormA + LB/quadFormB);

    vectorToProject = xJ + halfNablaF / halfL;

    % project onto the feasible set
    % Projection code by
    % Condat, Laurent. 
    % "Fast projection onto the simplex and the l1 Ball." 
    % Mathematical Programming 158.1 (2016): 575-585    
    uJ = max(vectorToProject-max((cumsum(sort(vectorToProject,2,'descend'),2)-sum(xJ)).*onesRowVectorInv),0);
  
    % compute dJ and the next iterate of x
    dJ = (uJ - xJ)';
    x(J) = uJ;

    % update the variables tracking x'*A*x and x'*B*x
    quadFormB = quadFormB + 2*(z_B*dJ)+ dJ'*subB*dJ; 
    quadFormA = quadFormA + 2*(z_A*dJ) +dJ'*subA*dJ;      

    iterCount=iterCount+1;
end

numIter = iterCount;

% output objective value
objValue = (x*A*x')/(x*B*x');

% unit test
if abs(sum(x) - 1) > 0.00001 || any(x < 0)
    error("Infeasible vector x!")
end

end
