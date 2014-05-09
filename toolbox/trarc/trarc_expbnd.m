function [stp bndflg] = trarc_expbnd(phi,x,a,bl,stpmax,infval)
  %trarc_expbnd  compute first intersection of bl <= x + a'*phi(stp)
  %
  % Computes the first point of intersection of
  %
  %  bl <= x + a'*phi(stp)
  %
  % with stp <= stpmax.
  %
  % It is assumed that bl < x + a'*phi(0) sufficiently so that a root is not
  % computed at stp = 0.
  %
  % Input:
  %  phi = trarc object (vector valued, length n)
  %  x = initial value (scalar)
  %  a = constraint vector (legnth n)
  %  bl = lower bound (scalar)
  %  stpmax = max step size
  %  infval = infinity value
  %
  % Output:
  %  stp = intersecting step size or stpmax
  %  bndflg = 1 if intersection exists
  %           0 if no intersection exists
  %
  
  % initialize output
  stp = stpmax;
  bndflg = 0;
  
  % check if bound is finite
  if bl > -infval
    
    % compute the intersection
    r = phi.intersect(a,bl-x);
    
    % check if intersection exists and update accordingly
    if ~isempty(r) && r <= stpmax
      stp = r;
      bndflg = 1;
    end
    
  end
  
end

