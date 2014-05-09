function [stp bndflg] = trarc_expstp(phi,x,a,bl,bu,stpmax,infval)
  %trarc_expstp  compute first intersection of bl <= x + a'*phi(stp) <= bu
  %
  % Computes the first point of intersection of
  %
  %  bl <= x + a'*phi(stp) <= bu
  %
  % with stp <= stpmax.
  %
  % It is assumed that bl < x + a'*phi(0) < bu sufficiently so that a
  % root is not computed at stp=0.
  %
  % Input:
  %  phi = trarc object (vector valued, length n)
  %  x = initial value (scalar)
  %  a = constraint vector (length n)
  %  bl = lower bound (scalar)
  %  bu = upper bound (scalar)
  %  stpmax = max step size
  %  infval = infinity value
  %
  % Output:
  %  stp = intersecting step size or stpmax
  %  bndflg = -1 if intersection at lower bound
  %            0 if no intersection exists
  %           +1 if intersection at upper bound
  %
  
  % compute intersections with lower and upper bounds
  [stp_bl bnd_bl] = trarc_expbnd(phi,x,a,bl,stpmax,infval);
  [stp_bu bnd_bu] = trarc_expbnd(phi,-x,-a,-bu,stpmax,infval);
  
  % return the smallest step size
  if bnd_bl == 0 && bnd_bu == 0
    % no bound is hit
    stp = stpmax;
    bndflg = 0;
  elseif stp_bl < stp_bu
    % lower bound is hit
    stp = stp_bl;
    bndflg = -1;
  else
    % upper bound is hit
    stp = stp_bu;
    bndflg = 1;
  end
  
end

