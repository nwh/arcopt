function [alpha r b] = expand_rt(idx,x,p,bl,bu,t,infval)
  %expand_rt  ratio test for computing largest possible step length
  %
  % Given variable vector x such that bl <= x <= bu and search
  % direction p, computes max step length alpha such that
  % bl <= x + alpha*p <= bu.
  %
  % If on completion, alpha = infval and j = 0, then the search direction
  % p is unbounded.
  %
  % Input:
  %  idx = indices to look at
  %  x = variable vector, length n
  %  p = search direction, length n
  %  bl = lower bound, length n
  %  bu = upper bound, length n
  %  t = tolerance designating negligible values of p
  %  infval = infinity value
  %
  % Output:
  %  alpha = max step length
  %  r = blocking variable
  %  b = -1 if l is hit
  %       0 if unbounded
  %       1 if u is hit
  %
  
  if nargin < 5
    t = 0.0;
  end
  
  alpha = infval;
  r = 0;
  b = 0;
  
  for i = idx(:)'
    [alpha_temp b_temp] = expand_stp(x(i),p(i),bl(i),bu(i),t,infval);
    if alpha_temp < alpha
      alpha = alpha_temp;
      b = b_temp;
      r = i;
    end
  end
  
end
