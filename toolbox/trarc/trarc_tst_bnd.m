%trarc_tst_bnd  test function for trarc's bound method
%
% The trarc.bound method returns steps sizes such that the solution norm has a
% particular value.  Say arc1 is the name of the trarc object.  If we call
%
%   s = arc1.bound(c);
%
% we should get norm(arc1.sol_s(s))^2 == c, where c > 0.  This can be used when
% it is desired that the trarc solution has a particular norm.  This function
% tests the trarc.bound method.
%
% Example:
%   n = 13;
%   H = randn(n,n);
%   H = 0.5*(H+H');
%   g = -randn(n,1);
%   c = 1.5;
%   [tst_flg wtw resid] = trarc_tst_bnd(H,g,c)
%
% Input:
%   H = Hessian matrix
%   g = gradient vector
%   c = input value to bound()
%   tol = error tolerance, defaults to 1e-10
%
% Output:
%   tst_flg = true if test passes, false otherwise
%   wtw = inner product of arc solution
%   resid = wtw - c
%

function [tst_flg wtw resid] = trarc_tst_bnd(H,g,c,tol)
  
  % set the optional tolerance
  if nargin < 6 || isempty(tol)
    tol = 1e-10;
  end
  
  % construct arc
  arc1 = trarc(H,g);

  % get the step size of hitting the bound
  s_bnd = arc1.bound(c);
  
  % evaluate arc at computed step size
  w_bnd = arc1.sol_s(s_bnd);

  % evaluate constraint and residual
  wtw = w_bnd'*w_bnd;
  resid = wtw - c;
  
  % set the tst_flg
  if abs(resid) <= tol
    tst_flg = true;
  else
    tst_flg = false;
  end
  
end