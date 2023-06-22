function M_tilde = computeM_EiC(A,B,xStar)
% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
%
% given solution vector x, this function checks whether x satisfies the KKT
% conditions. For the problem
% max x'Ax
% note that this is equivalent to
% min -x'Ax
% with gradient
% nabla f = -2Ax

% output is boolean, indicating whether or not x satisfies KKT conditions.
xStar = xStar(:);
n = numel(xStar);

linA = A*xStar; linB = B*xStar;
nablaF = -2*(linA / (xStar'*linA) - linB / (xStar'*linB) );

x = sdpvar(n,1);
obj = nablaF'*x - nablaF'*xStar;

cons = [x >= 0; sum(x) == 1];
ops = sdpsettings('verbose',0);
diagnostics = optimize(cons,obj,ops);
M_tilde = -value(obj);
end


