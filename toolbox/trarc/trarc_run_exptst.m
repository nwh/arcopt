%trarc_run_exptst  driver to test mgf_lc.arctest

function trarc_run_exptst
  
  % seed the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate arc data
  n = 10;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  
  % generate a constraint
  m = 20;
  A = randn(m,n);
  bl = -1-rand(m,1);
  bu = 1+rand(m,1);
  
  % generate starting point
  x = 2*rand(m,1)-1;
  
  % set delta
  delta = 1e-6;
  
  % max step size
  stpmax = .3;
  
  % tolerance
  tol = 1e-10;

  % construct arc
  arc = trarc(H,g);
  
  % run the test
  [tstflg stp bndflg] = trarc_tst_exptst(arc,x,A,bl,bu,stpmax,delta,tol);
  
  % print the results
  fprintf('\n... trarc exptst test ...\n')
  fprintf('Problem size = %d\n',n);
  fprintf('Number of constraints = %d\n',m);
  fprintf('bndflg = %g\n',bndflg);
  fprintf('stp = %g\n',stp);
  fprintf('stpmax = %g\n',stpmax);
  fprintf('test flag = %d\n',tstflg);
  
  if ~tstflg
    keyboard
  end
  
end