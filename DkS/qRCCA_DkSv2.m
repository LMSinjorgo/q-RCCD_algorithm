function [runningTime, objValue, x] = qRCCA_DkS(A,q,k,M,T,eps)
rng(123)
% Implementation of the qRCCA algorithm for 
% the continuous relaxation of DkS (densest k subgraph problem).
% Lennart Sinjorgo & Renata Sotirov

% DkS:
% minimize: x'Ax,
% subject to: sum(x) = k, 0 \leq x \leq 1. Here, A is the adjaency matrix
% of a graph

% Each step, select q coordinates (uniformly random) to update.

% For the QP:
% max f(u)  = NablaF (u - xJ) - L/2 * || u - xJ ||^2
%           = NablaF *u - L/2 * (u'u - 2u'x)      + Constants
%           = - u' *(L/2)* u + (NablaF + L*x)'*u  + Constants
%           = 0.5 u'*( -L*I )*u + c'u                  + Constants
% s.t. aJ * uJ = aJ * xJ

tic
n = size(A,2);

% initial feasible starting point
x = (k/n) * ones(n,1);
linA = A*x;
quadForm = x'*linA;

defaultopt = mskoptimset;
options = mskoptimset(defaultopt,[]);
[cmd,~,param] = msksetup(1,options);

prob.qosubi = (1:q)';
prob.qosubj = prob.qosubi;
prob.blx    = zeros(1,q);
prob.bux    = ones(1,q);
prob.a = sparse(ones(1,q));

onesVector = prob.bux';
objOld = quadForm;

m = 1;

% compute Lipschitz constant
L = full(max(sum(A,2)));
stagCounter = 1;
fprintf("  Iter \t obj Value \n");
fprintf("  0    \t %.3f \n",quadForm);
while m < M
    J = randperm(n,q);
    xJ = x(J);

    % compute partial derivatives
    NablaF = 2 * linA(J,:);
    
    % compute parameters
    sum_xJ = sum(xJ);
    linearCoeff = NablaF + L * xJ;
    
    % Set MOSEK parameters
    prob.qoval = L * onesVector;
    prob.c = -(linearCoeff);
    prob.blc = sum_xJ;
    prob.buc = sum_xJ;

    % solve the quadratic programme
    [~,res] = mosekopt(cmd,prob,param);
    uJ = res.sol.itr.xx;

    % compute dJ and the next iterate of x
    dJ = uJ - xJ;
    x(J) = uJ;
    
    % update the variables tracking Ax, x'Ax, Bx, x'Bx
    z_A = A(:,J)*dJ; 
    
    quadForm = quadForm + 2*dJ'*linA(J) + dJ'*z_A(J,:);
    linA = linA + z_A;
    
    % track the objective increase
    if (quadForm - objOld) < eps
        stagCounter = stagCounter+1;
        if stagCounter > T
            break;
        end
    else
        stagCounter = 0;
    end
    
    objOld = quadForm;
    
    m = m+1;
    fprintf("  %g    \t %.3f \n",m,quadForm);
end
runningTime = toc;

objValue = quadForm;
end
