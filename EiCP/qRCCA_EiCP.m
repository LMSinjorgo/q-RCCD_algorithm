function [runningTime, objValue, x] = qRCCA_EiCP(A,B,q,maxIter,k,eps)
% Implementation of the qRCCA algorithm for EiCP.
% Lennart Sinjorgo & Renata Sotirov

% EiCP:
% minimize: log(x'Ax) - log(x'Bx),
% subject to: sum(x) = 1, x \geq 0

% matrices A and B are sparse, symmetric, nonnegative and have positive
% diagonal

% Each step, select q coordinates (uniformly random) to update.

% For the QP:
% max f(u)  = NablaF (u - xJ) - L/2 * || u - xJ ||^2
%           = NablaF *u - L/2 * (u'u - 2u'x)      + Constants
%           = - u' *(L/2)* u + (NablaF + L*x)'*u  + Constants
%           = 0.5 u'*( -L*I )*u + c'u                  + Constants
% s.t. aJ * uJ = aJ * xJ

tic
n = size(A,2);

x = (1/n) * ones(n,1);
quadFormA = x'*A*x; quadFormB = x'*B*x;
linA = A*x; linB = B*x;


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


currentIter = 1;

objDiff = ones(k,1);

objDiffPlot = [];
while currentIter <  maxIter
    J = randperm(n,q);
    xJ = x(J);
    
    % compute Lipschitz constant
    L = 2*(quickNormEst(A(J,J),q)/quadFormA + quickNormEst(B(J,J),q)/quadFormB );    

    % compute partial derivatives
    NablaF = 2 * ( linA(J,:) / quadFormA - linB(J,:) / quadFormB );
    
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
    quadFormA = quadFormA + 2*dJ'*linA(J) + dJ'*z_A(J,:);
    quadFormB = quadFormB + 2*dJ'*linB(J) + dJ'*z_B(J,:);

    linA = linA + z_A; linB = linB + z_B;
    
    obj = log(quadFormA)-log(quadFormB);
    objDiff(mod(currentIter,k)+1) = obj - objOld;
    
    objDiffPlot = [objDiffPlot; obj - objOld];
    objOld = obj;
    
    
    currentIter = currentIter+1;
    
    if sum(objDiff) < eps
        break;
    end
end
runningTime = toc;

plot(objDiffPlot)
objValue = log(x'*A*x)-log(x'*B*x);
end

