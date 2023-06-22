% On convergence of a q-random coordinate constrained algorithm for non-convex problems (2023)
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
clear
clc

% set parameters
n = 3000;
d = 10^(-2);
q = 750;
maxTime = 10;
maxIter = Inf;

% generate matrices
A = generateMatrix(n,d);
B = generateMatrix(n,d);

% Provide no feasible starting point (q-RCCD algorithm will pick scaled all
% ones vector)
xStart = [];

% run q-RCCD algorithm
[objValue, x,iterCount] = qRCCD_EiC(A,B,q,maxTime,maxIter,xStart);

% check if x is a stationairy point
% (Requires YALMIP)
% (For M \approx 0, runTime should be approx 2/3 minutes)
M_tilde = computeM_EiC(A,B,x);

