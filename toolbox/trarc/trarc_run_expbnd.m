%trarc_run_expbnd  driver to test trarc_expbnd

function trarc_run_expbnd
  
  % seed the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate arc data
  n = 10;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  
  % generate a constraint
  a = randn(n,1);
  bl = -rand(1,1);
  
  % generate starting point
  x = rand(1,1);
  
  % max step size
  stpmax = 10;
  
  % tolerance
  tol = 1e-10;

  % construct arc
  arc = trarc(H,g);
  
  % run the test
  [tstflg resid xnew stp bndflg] = trarc_tst_expbnd(arc,x,a,bl,stpmax,tol);
  
  % print the results
  fprintf('\n... trarc arcbnd test ...\n')
  fprintf('Problem size = %d\n',n);
  fprintf('x = %g\n',x);
  fprintf('bl = %g\n',bl);
  fprintf('xnew = %g\n',xnew);
  fprintf('residual = %g\n',resid);
  fprintf('bndflg = %g\n',bndflg);
  fprintf('stp = %g\n',stp);
  fprintf('stpmax = %g\n',stpmax);
  fprintf('test flag = %d\n',tstflg);

  if ~tstflg
    keyboard
  end
  
end