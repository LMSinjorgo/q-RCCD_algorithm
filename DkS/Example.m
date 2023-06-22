% On convergence of a q-random coordinate constrained algorithm for non-convex problems
% By:
% A. Ghaffari-Hadigheh (Azarbaijan Shahid Madani University, Tabriz, Iran)
% L. Sinjorgo (EOR Department, Tilburg University, The Netherlands)
%           corresponding author: l.m.sinjorgo@tilburguniversity.edu
% R. Sotirov (EOR Department, Tilburg University, The Netherlands)
%                   year: 2023
clear
clc

% load adjacency matrix of graph
load("C1000-9.mat");

% or load another graph
%load("CA-AstroPh.mat");

% set parameters
q = 100;
maxTime = 10;
maxIter = Inf;
k = 200;

% Provide no feasible starting point (qRCCD algorithm will pick scaled all
% ones vector)
xStart = [];

% run q-RCCD algorithm
[objValue, x,iterCount] = qRCCD_DkS(A,q,k,maxTime,maxIter,xStart);

% extract feasible binary vector
[~,binaryIdx] = maxk(x,k);
subMat = A(binaryIdx,binaryIdx);
LB = sum(sum(subMat));

% visualization of the dense subgraph
plot(graph(subMat));

% check if x is a stationairy point
% Requires YALMIP (or reformulate with MATLAB linprog)
M_tilde = computeM_DkS(A,x,k);

