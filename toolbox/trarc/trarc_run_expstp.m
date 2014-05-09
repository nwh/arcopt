%trarc_run_expstp  driver to test trarc_expstp

function trarc_run_expstp
  
  % seed the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate arc data
  n = 10;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  
  % generate a constraint
  a = randn(n,1);
  bl = -1-rand(1,1);
  bu = 1+rand(1,1);
  
  % generate starting point
  x = 2*rand(1,1)-1;
  
  % max step size
  stpmax = 10;
  
  % tolerance
  tol = 1e-10;

  % construct arc
  arc = trarc(H,g);
  
  % run the test
  [tstflg resid_l resid_u xnew stp bndflg] = trarc_tst_expstp(arc,x,a,bl,bu,stpmax,tol);
  
  % print the results
  fprintf('\n... trarc expstp test ...\n')
  fprintf('Problem size = %d\n',n);
  fprintf('x = %g\n',x);
  fprintf('bl = %g\n',bl);
  fprintf('bu = %g\n',bu);
  fprintf('xnew = %g\n',xnew);
  fprintf('residual at bl = %g\n',resid_l);
  fprintf('residual at bu = %g\n',resid_u);
  fprintf('bndflg = %g\n',bndflg);
  fprintf('stp = %g\n',stp);
  fprintf('stpmax = %g\n',stpmax);
  fprintf('test flag = %d\n',tstflg);
  
  if ~tstflg
    keyboard
  end
  
end