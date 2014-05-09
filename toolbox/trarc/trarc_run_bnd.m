%trarc_run_bnd  basic testing script for trarc.bound
%
% This script will test the computation of the step length where the arc
% intersects the trust region bound.
%

function trarc_run_bnd
  
  % seed the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate arc data
  n = 10;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  
  % constraint value
  c = 2;
  
  % tolerance
  tol = 1e-10;
  
  % run the test
  [tst_flg wtw resid] = trarc_tst_bnd(H,g,c,tol);
  
  % print the results
  fprintf('\n... trarc bound test ...\n')
  fprintf('Problem size = %d\n',n);
  fprintf('c = %g\n',c);
  fprintf('w(s)''*w(s) = %g\n',wtw);
  fprintf('w(s)''*w(s) - c = %g\n',resid);
  fprintf('test flag = %d\n',tst_flg);
  
end