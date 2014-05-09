%trarc_run  basic testing script for trarc
%
% This script will test the evaluation of the arc and derivatives
%

function trarc_run
  
  % set the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate data
  n = 10;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  s = 3;
  
  % construct arc
  arc = trarc(H,g);
  
  % parameters
  h = 1e-8;
  htol = 1e-5;
  mtol = 1e-10;
  
  % run the test
  [tst_flg err_m err_fd1 err_fd2] = trarc_tst(H,g,s,h,htol,mtol);
  
  % print the results
  fprintf('\n... trarc basic test ...\n')
  fprintf('Problem size = %d\n',n);
  fprintf('max eigenvalue = %g\n',arc.v_max);
  fprintf('min eigenvalue = %g\n',arc.v_min);
  fprintf('w(s) error = %g\n',err_m);
  fprintf('w''(s) error = %g\n',err_fd1);
  fprintf('w''''(s) error = %g\n',err_fd2);
  fprintf('test flag = %d\n',tst_flg);
  
end