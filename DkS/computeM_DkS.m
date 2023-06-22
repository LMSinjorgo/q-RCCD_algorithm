function M_tilde = computeM_DkS(A,xStar,k)
% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
%
% REQUIRES YALMIP!
% given solution vector x, this function checks whether x satisfies the KKT
% conditions. For the problem
% max x'Ax
% note that this is equivalent to
% min -x'Ax
% with gradient
% nabla f = -2Ax

% output is boolean, indicating whether or not x satisfies KKT conditions.
n = numel(xStar);

nablaF = -2*A*xStar;

x = sdpvar(n,1);
obj = nablaF'*x - nablaF'*xStar;

cons = [x >= 0; x <= 1; sum(x) == k];
ops = sdpsettings('verbose',0);
diagnostics = optimize(cons,obj,ops);
M_tilde = -value(obj);
end


