function [runningTime, objValue, x] = qRCCA_EiCP(A,B,q,M,T,eps)
% Implementation of the qRCCA algorithm for EiCP.
% Lennart Sinjorgo & Renata Sotirov



% EiCP:
% minimize: log(x'Ax) - log(x'Bx),
% subject to: sum(x) = 1, x \geq 0

% inputs:
% matrices A and B are sparse, symmetric, nonnegative and have positive
% diagonal.
% Each step, select q coordinates (uniformly random) to update.
% M is maximum number of iterations
% STOP iterations in iteration m if (obj now) - (obj at step m-T) < eps

% For the QP:
% max f(u)  = NablaF (u - xJ) - L/2 * || u - xJ ||^2
%           = NablaF *u - L/2 * (u'u - 2u'x)      + Constants
%           = - u' *(L/2)* u + (NablaF + L*x)'*u  + Constants
%           = 0.5 u'*( -L*I )*u + c'u                  + Constants
% s.t. aJ * uJ = aJ * xJ

tic
n = size(A,2);

x = (1/n) * ones(n,1);
linA = A*x; linB = B*x;
quadFormA = x'*linA; quadFormB = x'*linB;

% Set up MOSEK Parameters
defaultopt = mskoptimset;
options = mskoptimset(defaultopt,[]);
[cmd,~,param] = msksetup(1,options);

prob.qosubi = (1:q)';
prob.qosubj = prob.qosubi;
prob.blx    = zeros(1,q);
prob.bux    = ones(1,q);
prob.a = sparse(ones(1,q));

onesVector = prob.bux';
objOld = log(quadFormA)-log(quadFormB);


m = 1;

% structures to check the stopping criteria
T = T+1; % matlab starts counting at 0.
objMatrix = -ones(1,T);
objMatrix(1) = objOld;

while m <  M
    J = randperm(n,q);
    xJ = x(J);
    
    % compute Lipschitz constant
    if q < 15
        % norm of full matrix is faster than normest for smaller matrices
        L = 2*(norm(full(A(J,J)))/quadFormA + norm(full(A(J,J))) / quadFormB );
    else
        L = 2*(quickNormEst(A(J,J),q)/quadFormA + quickNormEst(B(J,J),q)/quadFormB );    
    end

    % compute partial derivatives
    NablaF = 2 * ( linA(J) / quadFormA - linB(J) / quadFormB );
    
    % Set MOSEK parameters
    sum_xJ = sum(xJ);
    prob.blc    = sum_xJ;
    prob.buc    = sum_xJ;
    prob.qoval = L * onesVector;
    prob.c = -(NablaF + L * xJ);
    
    % solve the quadratic programme
    [~,res] = mosekopt(cmd,prob,param);
    uJ = res.sol.itr.xx;
    
    % compute dJ and the next iterate of x
    dJ = uJ - xJ;
    x(J) = uJ;
    
    % update the variables tracking Ax, x'Ax, Bx, x'Bx
    z_A = A(:,J)*dJ; 
    z_B = B(:,J)*dJ;
    quadFormA = quadFormA + 2*dJ'*linA(J) + dJ'*z_A(J);
    quadFormB = quadFormB + 2*dJ'*linB(J) + dJ'*z_B(J);

    linA = linA + z_A; linB = linB + z_B;
    
    % track the objective
    obj = log(quadFormA/quadFormB);

    modCounter1 = mod(m,T)+1;
    modCounter2 = mod(m+1,T)+1;
    objMatrix(1,modCounter1) = obj;

    if objMatrix(1,modCounter1)-objMatrix(1,modCounter2) < eps
        % objective has not increased with more than epsilon over the last
        % T iterations.
        break;
    end
    m = m+1;
end
runningTime = toc;

objValue = obj;
end
