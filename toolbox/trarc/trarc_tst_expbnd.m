%trarc_tst_expbnd test function for trarc_expbnd

function [tstflg resid xnew stp bndflg] = trarc_tst_expbnd(arc,x,a,bl,stpmax,tol)
  
  % set default tolerance if needed
  if nargin < 6 || isempty(tol)
    tol = 1e-10;
  end
  
  % set the infval
  infval = 1e20;
  
  % orient a
  a = a(:);
  
  % evaluate expbnd method
  [stp bndflg] = trarc_expbnd(arc,x,a,bl,stpmax,infval);
  
  % evaluate arc and compute residual
  w = arc.sol_s(stp);
  xnew = x+a'*w;
  resid = xnew - bl;

  % test the results
  if bndflg && abs(resid) <= tol
    tstflg = 1;
  elseif ~bndflg && resid > tol
    tstflg = 1;
  else
    tstflg = 0;
  end
    
end