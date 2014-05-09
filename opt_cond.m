function [optflag,opt_p,opt_dl,opt_du] = opt_cond(x,y,g,A,bl,bu,ptol,dtol)
%opt_cond  check first order optimality conditions for problem LC1
%
% This is a simple function to check the optimality condition of the
% optimization problem:
%
%   min f(x)
%   s/t A*x = s
%       bl <= [x; s] <= bu
%
% x are problem variables
% s are slack variables
%
% If we stack x and s such that xs = [x; s], the optimality conditions are:
% 
%   min(xs-bl,bu-xs) >= -ptol
%   min(xs-bl, z) <= dtol
%   min(bu-xs,-z) <= dtol
%
% ptol and dtol are numerical tolerances.
%
% Usage:
%   [optflag,opt_p,opt_dl,opt_du] = opt_cond(x,y,g,A,bl,bu,ptol,dtol)
%
% Input:
%   x = variable vector, length n (slack variables are computed)
%   y = multiplier vector, length m
%   g = function gradient, length n
%   A = constraint matrix, size m by n
%   bl = lower bound on all variables, length m+n
%   bu = upper bound on all variables, length m+n
%   ptol = primal tolerance (defaults to 1e-6)
%   dtol = dual tolerance (defaults to 1e-6)
%
% Output:
%   optflag = 1 if optimal, 0 if not optimal
%   opt_p = primal infeasibility check
%   opt_dl = lower bound dual infeasibility check
%   opt_du = upper bound dual infeasibility check
%

%
% 2010-11-15 first version.
%

if nargin < 7 || isempty(ptol)
  ptol = 1e-6;
end

if nargin < 8 || isempty(dtol)
  dtol = 1e-6;
end

% get problem size
n = length(x);
m = length(y);

% if there are no constraints, create dummy matrix
if m == 0
  A = zeros(0,n);
end

% compute reduced gradient
z = zeros(m+n,1);
if m > 0
  z(1:n) = g-A'*y;
  z(n+1:m+n) = y;
else
  z = g;
end
  
% compute constraint values
% these are also the value of the slack variables
s = A*x;

% augment to full variable vector
xs = [x; s];

% check feasibility (primal tolerance)
opt_p = min(min(xs-bl,bu-xs));

% check dual tolerances
opt_dl = max(min(xs-bl,z));
opt_du = max(min(bu-xs,-z));

if opt_p >= -ptol && opt_dl <= dtol && opt_du <= dtol
  % solution (x,y) is within first order optimality tolerance
  optflag = 1;
else
  % (x,y) not optimal
  optflag = 0;
end