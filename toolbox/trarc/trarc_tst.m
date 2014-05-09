%trarc_tst  test function for trust region arc function and derivatives
%
% This function provides some simple tests on the evaluation of the trust
% region arc and it's derivatives.
%
% The first test uses the arc class to compute the multiplier value:
%   y = 1/s - v_min(H)
% It then compares the solution of (H+y*I)*w = -g to the output from arc.sol_s.
% These values should be nearly equivalent.  If s = 0, it checks to make sure
% norm(arc.sol_s(0)) <= mtol.
%
% The second test uses a finite difference to check the analytic first
% derivative.
%
% The third test uses a finite difference to check the analytic second
% derivative.
% 
% Example:
%   n = 10;
%   H = randn(n,n);
%   H = 0.5*(H+H');
%   g = -randn(n,1);
%   s = 5;
%   [tst_flg err_m err_fd1 err_fd2] = trarc_tst(H,g,s)
%
% Input:
%   H = Hessian matrix
%   g = gradient vector
%   s = trial step size
%   h = finite difference parameter, default = 1e-6
%   htol = finite difference error tolerance, default = h*1000
%   mtol = multiplier error tolerance, default = 1e-10
%
% Output:
%   tst_flg = true if all tests passed, false otherwise
%   err_m = error for multiplier check
%   err_fd1 = error for finite difference check on first derivative
%   err_fd2 = error for finite difference check on second derivative
%

function [tst_flg err_m err_fd1 err_fd2] = trarc_tst(H,g,s,h,htol,mtol)
  
  % set optional tolerances
  if nargin < 4 || isempty(h)
    h = 1e-8;
  end
  
  if nargin < 5 || isempty(htol)
    htol = h*1000;
  end
  
  if nargin < 6 || isempty(mtol)
    mtol = 1e-10;
  end
  
  % get problem size
  n = length(g);
  
  % construct the arc
  arc1 = trarc(H,g);
  
  % test multiplier value
  % this test can only be done if s > 0.
  if s > 0
    y = arc1.mult(s);
    w_arc = arc1.sol_s(s);
    w_tr = -(H+y*eye(n))\g;
    err_m = norm(w_arc-w_tr);
  else
    w_arc = arc1.sol_s(s);
    err_m = norm(w_arc);
  end
  
  % test first derivative
  w = arc1.sol_s([s s+h]);
  w_td = arc1.sol_s(s,1);
  w_fd = (w(:,2) - w(:,1)) / h;
  err_fd1 = norm(w_fd-w_td);
  
  % test second derivative
  wd1 = arc1.sol_s([s s+h],1);
  w_td2 = arc1.sol_s(s,2);
  w_fd2 = (wd1(:,2) - wd1(:,1)) / h;
  err_fd2 = norm(w_fd2-w_td2);
  
  % set the output test flag
  if err_m <= mtol && err_fd1 <= htol && err_fd2 <= htol
    tst_flg = true;
  else
    tst_flg = false;
  end
  
end