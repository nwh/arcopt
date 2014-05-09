function [s b] = expand_stp(x,p,l,u,t,infval)
  %expand_stp  determine longest step for a single variable.
  %
  % This is a helper function to compute the longest step possible for
  % a single variable.  It guarantees that l <= x + s*p <= u.  Input t
  % is a tolerance designating negligible values of p.
  %
  % Input:
  %  x = current value
  %  p = search direction
  %  l = lower bound
  %  u = upper bound
  %  t = tolerance on p
  %  infval = infinity value
  %
  % Output:
  %  s = allowable step size
  %  b = -1 if l is hit
  %       0 if unbounded
  %       1 if u is hit
  %
  
  if nargin < 5
    t = 0.0;
  end
  
  if p < -t && l >= -infval
    s = (l-x)/p;
    b = -1;
  elseif p > t && u <= infval
    s = (u-x)/p;
    b = 1;
  else
    s = infval;
    b = 0;
  end
end
