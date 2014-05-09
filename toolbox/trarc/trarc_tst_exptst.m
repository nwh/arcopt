%trarc_tst_exptst  test function for trarc_exptst

function [tstflg stp bndflg] = ...
    trarc_tst_exptst(arc,x,A,bl,bu,stpmax,delta,tol)
  
  % set default tolerance if needed
  if nargin < 6 || isempty(tol)
    tol = 1e-10;
  end
  
  % set the infval
  infval = 1e20;
  
  % evaluate arcbound method
  [stp idx bndflg] = trarc_exptst(arc,x,A,bl,bu,stpmax,delta,infval);
    
  % evaluate arc and compute residual
  w = arc.sol_s(stp);
  xnew = x + A*w;
  resid_l = xnew - (bl - delta);
  resid_u = (bu + delta) - xnew;
  
  % test the results
  if bndflg == -1 && abs(resid_l(idx)) <= tol ... 
      && min(resid_l) >= -tol ... 
      && min(resid_u) >= -tol
    tstflg = 1;
  elseif bndflg == 1 && abs(resid_u(idx)) <= tol ... 
      && min(resid_l) >= -tol ... 
      && min(resid_u) >= -tol
    tstflg = 1;
  elseif ~bndflg && min(resid_l) >= -tol && min(resid_u) >= -tol
    tstflg = 1;
  else
    tstflg = 0;
  end

end