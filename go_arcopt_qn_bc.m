function go_arcopt_qn_bc
  
  % get options
  options = arcopt_qn_bc.optset();
  
  % algorithm options
  %options.search = 'line';
  %options.search = 'arc';
  %options.update = 'BFGS';
  %options.update = 'SR1';
  
  % step length selection options
  %options.ftol = 1e-5;
  %options.gtol = .9;
  %options.xtol = 1e-10;
  %options.fevmax_srch = 100;
  
  % general options
  %options.opttol = 1e-5;
  %options.fevmax = 1000;
  %options.itermax = 100;
  
  options.print_screen = 0;
  
  % test set cell
  test_set = {
    % name        n   bc setup           opttol
    'cvxqp_unc'   20  0  @cvxqp_setup    1e-6
    'cvxqp_con'   20  1  @cvxqp_setup    1e-6
    'rosen_unc'   2   0  @rosen_setup    1e-6
    'rosen_con'   2   1  @rosen_setup    1e-6
    'qp_con'      10  1  @qp_setup       1e-6
    'genrose_unc' 8   0  @genrose_setup  1e-6
    'genrose_con' 8   1  @genrose_setup  1e-6
    'var28_20'    20  1  @var28_setup    1e-4
    'var28_45'    45  1  @var28_setup    1e-4
    };
  
  [num_test junk] = size(test_set);
  
  % overall output
  qnbc_fid = fopen('output/qnbc_results.txt','w');
  
  h_str = '%12s %2s %2s %6s %8s %2s %4s %4s %9s %8s\n';
  d_str = '%12s %2d %2d %6g %8.2g %2d %4d %4d %9g %8d\n';
  fprintf(qnbc_fid,h_str,'name','n','bc','opttol','||Z''g||','ab','iter','nfev','wall_time','exitflag');
  
  test_index = 1:num_test;
  %test_index = 5;
  
  % run all tests
  for t = test_index
    name = test_set{t,1};
    n = test_set{t,2};
    bc = test_set{t,3};
    setup_func = test_set{t,4};
    prob_fid = fopen(['output/' name '.txt'],'w');
    opttol = test_set{t,5};
    
    [x0 lb ub n func] = setup_func(n);
    
    if bc == 0
      % problem is unconstrained
      % set very loose bounds
      lb = -100*ones(n,1);
      ub = 100*ones(n,1);
    end
    
    % set options
    options.ptol = opttol;
    options.print_file = prob_fid;
    
    fprintf(prob_fid,'Problem: %s\n',name);
    
    mysolver = arcopt_qn_bc(func,x0,lb,ub,options);
    [x info aux] = mysolver.solve();
    
    g = aux.g;
    ab = aux.xstate ~= 0;
    
    % collect output
    iter = aux.itercnt;
    nfev = aux.fevcnt;
    wall_time = aux.time;
    
    % print overall results
    fprintf(qnbc_fid,d_str,name,n,bc,opttol,norm(g(ab==0)),sum(abs(ab)),iter,nfev,wall_time,info);
    
    fclose(prob_fid);
  end
  fclose(qnbc_fid);
  
end
