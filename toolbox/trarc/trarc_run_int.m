%trarc_run_int  basic testing script for trarc.intersection
%
% This script will test the computation of the step length where the arc
% intersects a linear constraint.
%

function trarc_run_int
  
  % seed the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate arc data
  n = 10;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  
  % generate constraint data
  m = 20;
  A = randn(m,n);
  b = -rand(m,1);
  
  % max step size
  s_max = 1000;
  
  % tolerance
  tol = 1e-10;
  
  % run the test
  [tst_flg max_infeas s_int ix r_ix] = trarc_tst_int(H,g,A,b,s_max,tol);
  
  % print the results
  fprintf('\n... trarc intersection test ...\n')
  fprintf('Problem size = %d\n',n);
  fprintf('Number of constraints = %d\n',m);
  fprintf('max infeasibility = %g\n',max_infeas);
  fprintf('intersecting step size = %g\n',s_int);
  fprintf('index of limiting constraint = %d\n',ix);
  fprintf('residual with limiting constraint = %g\n',r_ix);
  fprintf('test flag = %d\n',tst_flg);
  
end