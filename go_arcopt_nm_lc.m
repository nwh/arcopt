% test rgd on cuter problem

function go_arcopt_nm_lc
  
  pdat_set = mcute_org2pdat('test-set.org');
  
  pset = [1 2 3 4 5 6 7];
  %pset = [1];
  
  pflg = zeros(size(pset));
  popt = zeros(size(pset));
  
  % test settings
  nmax = 300;
  h = 1e-8;
  tol = 1e-6;
  
  for p = pset
    
    %[pdat_test] = mcute_test_pdat(pdat_set(p),nmax,h,tol)
    %keyboard
    
    [prob build_flag build_msg] = mcute_build(pdat_set(p));
    [usrf usrh x0 bl bu A cl cu] = mcute_init(prob);
    
    %% set options
    options = arcopt_nm_lc.optset();
    options.fevmax = 10000;
    options.itermax = 2000;
    %options.crash = 'firstm';
    %options.facfrq = 0;
    %options.pcg_delta = 0;
    
    %% run mgf_lc
    mymgf_lc = arcopt_nm_lc(usrf,usrh,x0,bl,bu,A,cl,cu,options);
    [x f flg aux] = mymgf_lc.solve();
    
    pflg(p) = flg;
    
    %keyboard
    
    [optflag,opt_p,opt_dl,opt_du] = opt_cond(x,aux.y,aux.g,A,[bl;cl],[bu;cu]);
    
    popt(p) = optflag;
    
    %keyboard
    
  end
  
  fprintf('\n ... done solving ...\n');
  fprintf('prob flag opt_check\n');
  for p = pset
    fprintf('%4d %4d %d\n',p,pflg(p),popt(p));
  end
  
  % clean the directory
  mcute_clean;
  
end