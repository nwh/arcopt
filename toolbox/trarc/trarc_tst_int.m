%trarc_tst_int  test function for trarc's intersection method
%
% This function computes the intersection between an trust region arc and a set
% of linear inequality constraints.  It assumes that the initial point (w=0) is
% feasible.
%
% The function tests the final result for feasiblity.
%
% Example:
%   n = 3;
%   H = randn(n,n);
%   H = 0.5*(H+H');
%   g = -randn(n,1);
%   m = 3*n;
%   A = [eye(n); -eye(n); randn(n)];
%   b = -rand(m,1);
%   s_max = 100;
%   [tst_flg max_infeas s_int ix r_ix] = trarc_tst_int(H,g,A,b,s_max);
%
% Input:
%   H = Hessian
%   g = gradient
%   A = constraint matrix
%   b = constraint rhs
%   s_max = maximum step size
%   tol = infeasibility tolerance, defaults to 1e-10
%
% Output:
%   tst_flg = test flag, true if max_infeas <= tol
%   max_infeas = maximum infeasibility after evaluating at intersection point
%   s_int = step size for intersection, is s_max if no intersection
%   ix = index of limiting constraint, is 0 if no intersection
%   r_ix = residual of limiting constraint, is 0 if no intersection
%

function [tst_flg max_infeas s_int ix r_ix] = trarc_tst_int(H,g,A,b,s_max,tol)

  % set the optional tolerance
  if nargin < 6 || isempty(tol)
    tol = 1e-10;
  end
  
  % get the number of constraints
  m = length(b);
  
  % construct arc
  arc1 = trarc(H,g);
  
  % set initial values for step size and limiting constraint
  s_int = s_max;
  ix = 0;
  
  % evaluate intersection
  for i = 1:m
    s_int_i = arc1.intersect(A(i,:),b(i));
    
    if ~isempty(s_int_i) && s_int_i <= s_int
      s_int = s_int_i;
      ix = i;
    end
  end
  
  % compute resulting point
  w = arc1.sol_s(s_int);
  
  % check for infeasibility of resulting point
  resid = A*w-b;
  max_infeas = max(abs(resid(resid<0)));
  
  % get residual for limiting constraint
  if ix
    r_ix = resid(ix);
  else
    r_ix = 0;
  end
  
  % set final parameters
  if isempty(max_infeas)
    max_infeas = 0;
    tst_flg = true;
  elseif max_infeas <= tol
    tst_flg = true;
  else
    tst_flg = false;
  end
  
end