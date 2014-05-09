%trarc_tst_expstp  test function for trarc_expstp

function [tstflg resid_l resid_u xnew stp bndflg] = trarc_tst_expstp(arc,x,a,bl,bu,stpmax,tol)
  
  % set default tolerance if needed
  if nargin < 6 || isempty(tol)
    tol = 1e-10;
  end
  
  % set the infval
  infval = 1e20;
  
  % orient a
  a = a(:);
  
  % evaluate expstp method
  [stp bndflg] = trarc_expstp(arc,x,a,bl,bu,stpmax,infval);
  
  % evaluate arc and compute residual
  w = arc.sol_s(stp);
  xnew = x + a'*w;
  resid_l = xnew - bl;
  resid_u = bu - xnew;
  
  % test the results
  if bndflg == -1 && abs(resid_l) <= tol && resid_u >= -tol
    tstflg = 1;
  elseif bndflg == 1 && abs(resid_u) <= tol && resid_l >= -tol
    tstflg = 1;
  elseif ~bndflg && resid_l >= -tol && resid_u >= -tol
    tstflg = 1;
  else
    tstflg = 0;
  end

end