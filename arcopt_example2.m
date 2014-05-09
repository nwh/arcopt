% run single problem with arcopt_qn_bc

function arcopt_example2
  
  % setup problem
  [x0 xlow xupp n func] = genrose_setup(20);
  
  % solver options
  options = arcopt_qn_bc.optset();
  options.ftol = 1e-4;
  options.ptol = 1e-6;
  options.search = 'arc';
  options.update = 'SR1';
  %options.itermax = 15;
  
  % open debug file
  options.print_dbg = fopen('debug.org','w');
  %options.print_dbg = 1;
  
  % run solver
  mysolver = arcopt_qn_bc(func,x0,xlow,xupp,options);
  %mysolver = arcopt_qn_bc(func,x0,[],[],options);
  [xstar info aux] = mysolver.solve();
  
  % close debug file
  fclose(options.print_dbg);
  
  %keyboard
  
end
