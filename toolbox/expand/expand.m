function [alpha r b flg] = expand(idx,x,p,bl,bu,t,delta,tau,infval)
  %expand  expand procedure from Gill, Murray, and Saunders
  %
  % Computes step length alpha and blocking variable r subject to
  % parameters t, delta, and tau.
  %
  % See P. E. Gill et al., A Practical Anti-Cycling Procedure for
  % Linearly Constrained Optimization.  Mathematical Programming 45
  % (1989) 437-474
  %
  % Input:
  %  idx = indices to look at
  %  x = variable vector, length n
  %  p = search direction, length n
  %  bl = lower bound, length n
  %  bu = upper bound, length n
  %  t = tolerance designating negligible values of p
  %  delta = feasibility expansion
  %  tau = min step parameter
  %  infval = infinity value
  %
  % Output:
  %  alpha = step length
  %  r = blocking variable
  %  b = -1 if l is hit
  %       0 if unbounded
  %       1 if u is hit
  %  flg = 0 if normal step length is taken
  %      = 1 if computed step length is too small, alpha_min is taken
  %
  
  flg = 0;
  
  [alpha1 r1 b1] = expand_rt(idx,x,p,bl-delta,bu+delta,t,infval);
  r = 0;
  b = 0;
  pmax = 0;
  for j = idx(:)'
    [alpha b_temp] = expand_stp(x(j),p(j),bl(j),bu(j),t,infval);
    if alpha <= alpha1 && abs(p(j)) > pmax
      r = j;
      b = b_temp;
      alpha2 = alpha;
      pmax = abs(p(j));
    end
  end
  
  if r ~= 0
    % step is bounded, now check for degeneracy
    alpha_min = tau / abs(p(r));
    alpha = max(alpha2,alpha_min);
    
    if alpha2 < alpha_min
      % degenerate step
      flg = 1;
    end
  end
  
end
