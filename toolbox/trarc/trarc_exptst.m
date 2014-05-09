function [stp idx bndflg] = trarc_exptst(phi,x,A,bl,bu,stpmax,delta,infval)
  %trarc_exptst  compute first intersection of arc and linear constraints
  %
  % Computes the first point of intersection of
  %
  %  bl-delta <= x + A*phi(stp) <= bu+delta
  %
  % It is assumed that bl-delta < x + A*phi(0) < bu+delta sufficiently
  % so that a root is not computed at stp=0.
  %
  % Input:
  %  phi = trarc object (vector valued)
  %  x = initial value (vector)
  %  A = matrix
  %  bl = lower bound (vector)
  %  bu = upper bound (vector)
  %  stpmax = max step size
  %  delta = bound purturbation
  %  infval = infinity value
  %
  % Output:
  %  stp = intersecting step size or stpmax
  %  idx = index of limiting bound (0 if no bound hit)
  %  bndflg = -1 if intersection at lower bound
  %            0 if no intersection exists
  %           +1 if intersection at upper bound
  %
  %
  
  % set initial values
  stp = stpmax;
  idx = 0;
  bndflg = 0;
  
  % get size of x
  n = length(x);
  
  % loop through all variables testing for intersection with bound
  for i = 1:n
    
    % see if row i is limiting
    [stp_i bndflg_i] = trarc_expstp(phi,x(i),A(i,:),bl(i)-delta,bu(i)+delta,stpmax,infval);
    
    % if bounded, update flags
    if stp_i <= stp && bndflg_i ~= 0
      stp = stp_i;
      idx = i;
      bndflg = bndflg_i;
    end
    
  end
  
end

