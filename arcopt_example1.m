%mgf_lc_example  an indefinite QP as example for mgf_lc

function arcopt_example1
  
  % seed the rng
  rng_seed = 10;
  rng(rng_seed);

  % choose problem size
  n = 50; % number of variables
  m = 10; % number of constraints
  dQ = .1; % density for hessian
  dA = .1; % density for constraints
  
  % generate matrices
  Q = sprandsym(n,dQ);
  A = sprandn(m,n,dA);
  
  % generate special g
  [V D] = eig(full(Q));
  evals = diag(D);
  c = randn(n,1);
  c = V(:,evals >= 0)*c(evals >= 0);
  
  % generate rhs b
  x1 = rand(n,1);
  b = A*x1;
  
  % generate bounds
  bl = zeros(n,1);
  bu = ones(n,1);
  
  % choose starting point
  x0 = 0.5*ones(n,1);
  
  % get function handles
  usrfun = @(x) obj_func(x,Q,c);
  usrhess = @(x) obj_hess(x,Q);
  
  % get mgf_lc options structure
  mgf_opts = arcopt_nm_lc.optset();
  
  % set up solver object
  my_solver = arcopt_nm_lc(usrfun,usrhess,x0,bl,bu,A,b,b,mgf_opts);
  [x f flg aux] = my_solver.solve();
  
  %keyboard
  
end

function [f g] = obj_func(x,Q,c)
  
  f = 0.5*x'*(Q*x) + c'*x;
  g = Q*x + c;
  
end

function H = obj_hess(x,Q)
  
  H = Q;
  
end