classdef arcopt_qn_bc < handle
  
  properties
    opt = []; % solver options
    n = []; % number of variables in problem
    x0 = []; % initial point
    x = []; % variable vector
    xlow = []; % lower bound
    xupp = []; % upper bound
    usrfcn = []; % function handle
    
    f = []; % function value
    g = []; % gradient
    rg = []; % residual/reduced gradient vector
    f_old = []; % old function value
    g_old = []; % old gradient
    x_old = []; % old variable
    
    xstate = []; % variable state vector
    % xstate(i) = 2 if x(i) is fixed at xupp(i) for an iteration
    % xstate(i) = -1 if lower bound for x(i) is active
    % xstate(i) = 0 if x(i) is free
    % xstate(i) = 1 if upper bound for x(i) is active
    % xstate(i) = 2 if x(i) is fixed at xupp(i) for an iteration
    % xstate(i) = 3 if x(i) is fixed because xlow(i) == xupp(i)
    
    xfree = []; % logical vector indicating free variables
    xhit = []; % logical index indicating if variable hit bound
    
    B = []; % quasi-newton matrix
    p = []; % search direction
    
    stp = []; % step size
    stpmax = []; % max step size for iteration
    stpinit = []; % initial step size for iteration
    stpidx = []; % index of limiting variable in iteration
    stpbnd = []; % -1 if xlow hit, 1 if xupp hit, 0 if no bound hit
    
    info = []; % termination information
    info_srch = []; % line/arc search information
    
    time_start = []; % time at start
    time = []; % total wall time
    
    %% arc properties
    arc = []; % arc object
    arc_dim = []; % arc dimension
    arc_Q = []; % arc subspace matrix
    arc_H = []; % arc matrix term
    arc_g = []; % arc inhomogeneous term
    
    %% expand variables
    expcnt = 0; % expand counter
    exptol = 0; % working feasibility tolerance
    expinc = 0; % expand increase value (tau in paper)

    %% degeneracy
    dgnflg = 0; % degenerate step flag
    dgncnt = 0; % degenerate step counter
    
    %% counters
    fevcnt = []; % function evaluation counter
    itercnt = []; % iteration counter
    prntcnt = []; % print counter
    
  end
  
  methods (Static)
    
    function options = optset(varargin)
      % process options
      in_parse = inputParser;

      % general options
      in_parse.addParamValue('ptol',1e-6,@(x) x>0);
      in_parse.addParamValue('fevmax',10000,@(x) x>0);
      in_parse.addParamValue('itermax',1000,@(x) x>0);
      in_parse.addParamValue('infval',1e20,@(x) x>0);
      in_parse.addParamValue('init_scale',5,@(x) x>0);
      
      % aglorithm options
      in_parse.addParamValue('search','line',@(x) ismember(x,{'line','arc'}));
      in_parse.addParamValue('update','BFGS',@(x) ismember(x,{'BFGS','SR1'}));
      
      % lu fac options
      in_parse.addParamValue('lucond_tol',1e12,@(x) x>0);
      in_parse.addParamValue('lucond_mod',1e-6,@(x) x>0);
      
      % arc options
      in_parse.addParamValue('arc_vtol',1e-10,@(x) x > 0); % parameter for arc initial step size selection
      
      % search options
      in_parse.addParamValue('ftol',1e-4,@(x) x>0 && x<1); % arc search descent parameter
      in_parse.addParamValue('gtol',.9,@(x) x>0 && x<1); % arc search curvature parameter
      in_parse.addParamValue('xtol',1e-8,@(x) x>0); % arc search interval tolerance
      in_parse.addParamValue('fevmax_srch',100,@(x) x>0); % arc search max function evaluations
      in_parse.addParamValue('stpmax',1e10,@(x) x>0); % arc search max step size
      in_parse.addParamValue('stpmin',0,@(x) x>=0); % arc search min step size
      in_parse.addParamValue('bisect_srch',0,@(x) x==0 || x==1); % choose bisection for arc search
      
      % SR1 update options
      in_parse.addParamValue('sr1_tol',1e-8,@(x) x>0 && x<1);
      
      % print options
      in_parse.addParamValue('print_level','iter',@(x) ismember(x,{'iter','stats','none'}));
      in_parse.addParamValue('print_screen',1,@(x) ismember(x,[0 1]));
      in_parse.addParamValue('print_file',0,@(x) x == 0 | x > 2);
      in_parse.addParamValue('print_freq',10,@(x) x>0);
      in_parse.addParamValue('print_dbg',0,@(x) x >= 0);
      
      % expand parameters
      in_parse.addParamValue('expfrq',10000,@(x) x>=0); % expand reset frequency
      in_parse.addParamValue('expa',.5,@(x) x>0 && x<1); % parameter for initial feasibility tolerance
      in_parse.addParamValue('expb',.99,@(x) x>0 && x<1); % parameter for final feasibility tolerance
      in_parse.addParamValue('expsml',1e-11,@(x) x>=0); % tolerance for near zero numbers
      
      in_parse.parse(varargin{:});
      options = in_parse.Results;
    end
    
  end
  
  methods
    
    %% constructor
    function obj = arcopt_qn_bc(usrfcn,x0,xlow,xupp,varargin)
      
      % process options
      obj.opt = obj.optset(varargin{:});
      
      % process inputs
      if ~isvector(x0)
        error('arcopt:input','x0 is not a vector.');
      end
      
      % get number of variables
      obj.x0 = x0(:);
      obj.x = obj.x0;
      obj.n = length(obj.x);
      
      % process xlow
      if isempty(xlow)
        obj.xlow = -obj.opt.infval*ones(obj.n,1);
      elseif isscalar(xlow)
        obj.xlow = xlow*ones(obj.n,1);
      elseif ~isvector(xlow)
        error('arcopt:input','xlow must be empty, scalar or vector.')
      else
        if length(xlow) ~= obj.n
          error('arcopt:input','xlow has wrong size.');
        end
        obj.xlow = xlow(:);
      end

      % process xupp
      if isempty(xupp)
        obj.xupp = obj.opt.infval*ones(obj.n,1);
      elseif isscalar(xupp)
        obj.xupp = xupp*ones(obj.n,1);
      elseif ~isvector(xupp)
        error('arcopt:input','xupp must be empty, scalar or vector.')
      else
        if length(xupp) ~= obj.n
          error('arcopt:input','xupp has wrong size.');
        end
        obj.xupp = xupp(:);
      end

      % get function handle
      if ~isa(usrfcn,'function_handle')
        error('arcopt:input','usrfcn is not a function handle.')
      end
      
      % store function handle
      obj.usrfcn = usrfcn;
      
    end
    
    %% print methods
    function mprint(obj,level,varargin)
      %mprint  print handler
      %
      % This method provides printing functionality for the algoithms.
      %
      % Print level 'iter' is intended for printing each iteration
      % Print level 'stats' is intended for printing overall stats before
      % and after the algorithm runs.
      %
      % Usage:
      %   mprint('iter',[format string],[data])
      %   mprint('stats',[format string],[data])
      %
      % On entry:
      %   opt.print_level = print level desired by the user
      %   opt.print_screen = if true, mprint will print to the console
      %   opt.print_file = if greater than 0, is file id for printing
      %
      
      if ~ismember(level,{'iter','stats'})
        error('arcopt:mprint','invalid print level string: %s.',level);
      end
      
      % do nothing if print_level is 'none'
      if strcmp(obj.opt.print_level,'none')
        return;
      end
      
      % do nothing if print_level is 'stats' and level is 'iter'
      if strcmp(obj.opt.print_level,'stats') && strcmp(level,'iter')
        return;
      end
      
      % print to screen if requested
      if obj.opt.print_screen
        fprintf(1,varargin{:});
      end
      
      % print to file if requested
      if obj.opt.print_file
        fprintf(obj.opt.print_file,varargin{:});
      end
      
    end
    
    function dprint(obj,routine,message,value_fmt,value)
      % print to file if requested
      if obj.opt.print_dbg
        str_fmt = sprintf('| %%d | %%s | %%s | %s |\n',value_fmt);
        fprintf(obj.opt.print_dbg,str_fmt,obj.itercnt,routine,message,value);
      end
    end
    
    function print_start(obj)
      
      % is x0 in bounds
      if any(obj.x0 ~= obj.x)
        x0_test = 'out of bounds';
      else
        x0_test = 'in bounds';
      end
      
      obj.mprint('stats','qnbc problem information:\n')
      obj.mprint('stats','  number of variables: %d\n',obj.n)
      obj.mprint('stats','  number of fixed variables: %d\n',sum(obj.xstate == 3))
      
      obj.mprint('stats','\n')
      obj.mprint('stats','initial point information:\n')
      obj.mprint('stats','  number of initial free variables: %d\n',sum(obj.xfree))
      obj.mprint('stats','  x0: %s\n',x0_test)
      
      obj.mprint('stats','  max(|x|)....... = %g\n',max(abs(obj.x)))
      obj.mprint('stats','  f(x)........... = %g\n',obj.f)
      obj.mprint('stats','  max(|g(x)|).... = %g\n',max(abs(obj.g)))
      obj.mprint('stats','  max(|g(xfree)|) = %g\n',max(abs(obj.g(obj.xfree))))
      
    end
    
    function print_iter(obj)
      
      if ~strcmp(obj.opt.print_level,'iter')
        % do nothing
        return;
      end
      
      % currently printing
      % item: iter fevcnt f     ||rg|| stpinit  stp   d   srch |xfree|
      % head: iter fevcnt f     ||rg|| stpinit  stp   d   srch |xfree|
      % hspc: %5s  %6s    %8s   %7s    %7s      %7s   %1s %4s  %7s
      % dspc: %5d  %6d    %8.1e %7.1e  %7.1e    %7.1e %1s %4d  %7d
      
      if isempty(obj.prntcnt) || obj.prntcnt == obj.opt.print_freq
        obj.prntcnt = 1;
      else
        obj.prntcnt = obj.prntcnt + 1;
      end
      
      if obj.prntcnt == 1
        % print header
        obj.mprint('iter','\n');
        obj.mprint('iter','%5s  %6s  %8s  %7s  %7s  %7s  %1s  %4s  %7s\n',...
          'iter','fevcnt','f','||rg||','stpinit','stp','d','srch','|xfree|');
      end
      
      if obj.dgnflg
        dgnstr = 'd';
      else
        dgnstr = '';
      end
      
      obj.mprint('iter','%5d  %6d  %8.1e  %7.1e  %7.1e  %7.1e  %1s  %4d  %7d\n',...
        obj.itercnt,...
        obj.fevcnt,...
        obj.f,...
        max(abs(obj.g(obj.xfree))),...
        obj.stpinit, ...
        obj.stp,...
        dgnstr, ...
        obj.info_srch,...
        sum(obj.xfree));
      
    end
    
    function print_final(obj)
      
      obj.mprint('stats','\noptimization terminated:');
      
      obj.mprint('stats','\n  info: %d\n',obj.info);
      switch obj.info
        case 1
          obj.mprint('stats','  optimality conditions satisfied\n');
        case 2
          obj.mprint('stats','  reached maximum number of function evaluations\n');
        case 3
          obj.mprint('stats','  reached maximum number of iterations\n');
        case 4
          obj.mprint('stats','  search failure\n');
        otherwise
          obj.mprint('stats','  unknown termination condition\n');
      end
      
      obj.mprint('stats','\noptimization summary:\n');
      obj.mprint('stats','  number of iterations: %d\n',obj.itercnt)
      obj.mprint('stats','  number of function evaluations: %d\n',obj.fevcnt)
      obj.mprint('stats','  max(|x|)....... = %g\n',max(abs(obj.x)))
      obj.mprint('stats','  f(x)........... = %g\n',obj.f)
      obj.mprint('stats','  max(|g(x)|).... = %g\n',max(abs(obj.g)))
      obj.mprint('stats','  max(|g(xfree)|) = %g\n',max(abs(obj.g(obj.xfree))))
      obj.mprint('stats','  wall time (s).. = %g\n',obj.time)
      
      
    end
    
    %% main solve call
    function [xstar info aux] = solve(obj)
      
      % check method consistency
      if strcmp(obj.opt.search,'line')
        obj.opt.update = 'BFGS';
      end
      
      % project x if needed
      obj.x = max(obj.x,obj.xlow);
      obj.x = min(obj.x,obj.xupp);
      
      % initial function evaluation
      [obj.f obj.g] = obj.usrfcn(obj.x);
      obj.fevcnt = 1;
      
      % store old values
      obj.f_old = obj.f;
      obj.g_old = obj.g;
      obj.x_old = obj.x;
      
      % intialize xstatus
      obj.xstate = zeros(obj.n,1);
      
      % set up xstate
      obj.xstate( obj.x==obj.xlow & obj.g>0 ) = -1;
      obj.xstate( obj.x==obj.xupp & obj.g<0 ) = 1;
      obj.xstate( obj.xlow == obj.xupp ) = 3;
      
      % initialize xfree
      obj.xfree = obj.xstate == 0;
      
      % initialize qn matrix
      obj.B = (norm(obj.g)/obj.opt.init_scale)*eye(obj.n);
      
      % initialize counters
      obj.itercnt = 1;
      obj.dgncnt = 0;
      
      % initialize info
      obj.info = 0;
      
      % intialize expand
      obj.exp_init();
      
      % initial print
      obj.print_start();
      
      % start the timer
      obj.time_start = tic;
      
      % main optimization loop
      while(1)
        
        %if obj.itercnt >= 16
        %  keyboard
        %end
        
        %if obj.x(18) == obj.xupp(18)
        %  keyboard
        %end
        
        % check termination conditions
        obj.term_check();
        if obj.info
          break;
        end
        
        % run expand
        obj.exp_main();
        
        % compute reduced/residual gradient
        obj.comp_rg();
        
        % compute search direction
        obj.comp_p();
        
        % compute max step size
        obj.comp_stpmax_line();
        
        % compute arc
        obj.comp_arc();
        
        % compute stpmax along the arc
        obj.comp_stpmax_arc();
        
        % compute initial step size
        obj.comp_stpinit();
        
        % line search
        obj.comp_stp();
        
        % update B
        obj.update_B();
        
        % update state information
        obj.update_state();
        
        % some debug stuff
        if obj.opt.print_dbg
          obj.dprint('solve','norm(x-x_old)','%8.2e',norm(obj.x-obj.x_old,inf));
          obj.dprint('solve','norm(g-g_old)','%8.2e',norm(obj.g-obj.g_old,inf));
          obj.dprint('solve','f_old-f','%9.2e',obj.f_old-obj.f);
          if obj.f > obj.f_old
            obj.dprint('solve','exception','%s','f>f_old');
          end
        end
        
        % save information
        obj.x_old = obj.x;
        obj.f_old = obj.f;
        obj.g_old = obj.g;
        
        % print iteration information
        obj.print_iter();
        
        % update iteration counter
        obj.itercnt = obj.itercnt + 1;
        
        %obj
        %keyboard
        
      end
      
      % stop the timer
      obj.time = toc(obj.time_start);
      
      % print final information
      obj.print_final();
      
      % handle output
      xstar = obj.x;
      info = obj.info;
      aux.f = obj.f;
      aux.g = obj.g;
      aux.xstate = obj.xstate;
      aux.itercnt = obj.itercnt;
      aux.fevcnt = obj.fevcnt;
      aux.dgncnt = obj.dgncnt;
      aux.info = obj.info;
      aux.info_srch = obj.info_srch;
      aux.time = obj.time;
      
      %keyboard
      
    end
    
    %% helper methods
    function term_check(obj)
      
      if obj.fevcnt >= obj.opt.fevmax
        obj.info = 2;
      end
      
      if obj.itercnt >= obj.opt.itermax
        obj.info = 3;
      end
      
      if ismember(obj.info_srch,[2 3 4 7])
        obj.info = 4;
      end
      
      term_val = norm(obj.g(obj.xfree),'inf');
      term_tol = obj.opt.ptol*max(1,norm(obj.x,'inf'));
      if term_val <= term_tol && sum(obj.xhit) == 0
        obj.info = 1;
      end
      
      if obj.opt.print_dbg
        obj.dprint('term_check','term_val','%8.2e',term_val)
        obj.dprint('term_check','term_tol','%8.2e',term_tol)
        obj.dprint('term_check','sum(obj.xhit)','%d',sum(obj.xhit))
        obj.dprint('term_check','info','%d',obj.info)
      end
      
    end
    
    function comp_rg(obj)
      obj.rg = zeros(obj.n,1);
      obj.rg(obj.xfree) = obj.g(obj.xfree);
      
      if obj.opt.print_dbg
        obj.dprint('comp_rg','norm(g,inf)','%8.2e',norm(obj.g,inf));
        obj.dprint('comp_rg','norm(rg,inf)','%8.2e',norm(obj.g(obj.xfree),inf));
      end
      
    end
    
    function comp_p(obj)

      [L U pidx] = obj.comp_LU();
      
      obj.p = zeros(obj.n,1);
      gf = obj.g(obj.xfree);
      
      obj.p(obj.xfree) = -U\(L\gf(pidx));
      
      if obj.opt.print_dbg
        obj.dprint('comp_p','norm(p,inf)','%8.2e',norm(obj.p,inf));
        obj.dprint('comp_p','gtp','%9.2e',obj.g'*obj.p);
        obj.dprint('comp_p','gtp_normalized','%9.2e',obj.g'*obj.p/(norm(obj.g)*norm(obj.p)));
      end
      
    end
    
    function [L U pidx] = comp_LU(obj)
      %comp_LU  compute LU decomp of B
      
      lu_flg = 0;
      lucnt = 0;
      
      while lu_flg == 0
      
        % compute the lu decomp
        [L U pidx] = lu(obj.B(obj.xfree,obj.xfree),'vector');
        lucnt = lucnt + 1;
        
        ud = diag(U);
        umax = max(abs(ud));
        umin = min(abs(ud));
        
        lucond = umax / umin;
      
        if lucond >= obj.opt.lucond_tol
          obj.B = obj.B + eye(obj.n)*max(obj.opt.lucond_mod,obj.opt.lucond_mod*umax);
        else
          lu_flg = 1;
        end
      
      end

      if obj.opt.print_dbg
        obj.dprint('comp_LU','lucond','%9.2e',lucond);
        obj.dprint('comp_LU','umax','%9.2e',umax);
        obj.dprint('comp_LU','umin','%9.2e',umin);
        obj.dprint('comp_LU','lucnt','%9.2e',lucnt);
      end

    end
    
    function comp_arc(obj)

      % do nothing if method is line search
      if strcmp(obj.opt.search,'line')
        return;
      end
      
      % construct the arc

      % compute number of free variables
      nf = sum(obj.xfree);

      % construct the subspace matrix
      if nf == 1
        % only one free variable
        % use reduced gradient as subspace
        obj.arc_dim = 1;
        D = obj.rg;
      elseif nf >= 2
        % two free variables
        %use reduced gradient and dnc as subspace
        obj.arc_dim = 2;
        D = [obj.rg obj.p];
      end

      % compute orthogonal subspace matrix
      [obj.arc_Q R] = qr(D,0);

      % construct arc matrix
      obj.arc_H = obj.arc_Q'*(obj.B*obj.arc_Q);

      % construct arc rhs
      obj.arc_g = obj.arc_Q'*obj.rg;

      % create arc
      % (TODO) reuse same arc object
      obj.arc = trarc(obj.arc_H,obj.arc_g);
      
      if obj.opt.print_dbg
        obj.dprint('comp_arc','arc_dim','%d',obj.arc_dim);
        obj.dprint('comp_arc','min_eigval','%9.2e',obj.arc.v_min);
        obj.dprint('comp_arc','max_eigval','%9.2e',obj.arc.v_max);
      end
      
    end
    
    function comp_stpmax_line(obj)
 
      switch obj.opt.search
        case 'line'
          dx = obj.p;
        case 'arc'
          dx = -obj.rg;
      end
      
      idx = find(dx ~= 0);
      [obj.stpmax obj.stpidx obj.stpbnd obj.dgnflg] = expand(idx,...
        obj.x, ...
        dx, ...
        obj.xlow, ...
        obj.xupp, ...
        obj.opt.expsml, ...
        obj.exptol, ...
        obj.expinc, ...
        obj.opt.infval);

      % must check against stpmax parameter
      if obj.stpmax > obj.opt.stpmax
        % stpmax is too large, a bound will not be hit
        obj.stpmax = obj.opt.stpmax;
        obj.stpbnd = 0;
        obj.stpidx = 0;
        obj.dgnflg = 0;
      end
      
      % increment degeneracy counter if needed
      if obj.dgnflg
        obj.dgncnt = obj.dgncnt + 1;
      end
      
      % print debug info
      if obj.opt.print_dbg
        obj.dprint('comp_stpmax_line','stpmax','%9.2e',obj.stpmax);
        obj.dprint('comp_stpmax_line','stpidx','%d',obj.stpidx);
        obj.dprint('comp_stpmax_line','stpbnd','%d',obj.stpbnd);
        obj.dprint('comp_stpmax_line','dgnflg','%d',obj.dgnflg);
      end
      
    end
    
    function comp_stpmax_arc(obj)
      
      % do nothing if we are doing line search or if we are at a degenerate
      % point
      if strcmp(obj.opt.search,'line') || obj.dgnflg
        return;
      end
      
      [obj.stpmax obj.stpidx obj.stpbnd] = trarc_exptst(obj.arc, ...
        obj.x, ...
        obj.arc_Q, ...
        obj.xlow, ...
        obj.xupp, ...
        obj.opt.stpmax, ...
        obj.exptol, ...
        obj.opt.infval);
      
      if obj.opt.print_dbg
        obj.dprint('comp_stpmax_arc','stpmax','%9.2e',obj.stpmax);
        obj.dprint('comp_stpmax_arc','stpidx','%d',obj.stpidx);
        obj.dprint('comp_stpmax_arc','stpbnd','%d',obj.stpbnd);
      end
      
    end
    
    function comp_stpinit(obj)
      
      switch obj.opt.search
        case 'line'
          obj.stpinit = 1;
        case 'arc'
          if obj.arc.v_min >= obj.opt.arc_vtol
            % arc is in positive definite subspace
            obj.stpinit = 1/obj.arc.v_min;
          elseif obj.arc.v_min <= -obj.opt.arc_vtol
            % arc is in an indefinite subspace
            obj.stpinit = -1/obj.arc.v_min;
          else
            % arc is in positive semi-definite subspace
            obj.stpinit = 1;
          end
      end

      obj.stpinit = min(obj.stpinit,obj.stpmax);
      
      if obj.opt.print_dbg
        obj.dprint('comp_stpinit','stpinit','%9.2e',obj.stpinit);
      end
      
    end
    
    function [s ds] = arc_eval(obj,stp)
      %arc_eval parameterized arc evaluation for use with arc search
      %
      % Input:
      %   stp = trial step size
      %
      % On entry:
      %   arc = trust region arc object
      %   arc_Q = arc projection matrix
      %
      % Output:
      %   s = arc position at stp
      %   d = derivative of arc at stp
      %

      s = obj.arc_Q*obj.arc.sol_s(stp);
      ds = obj.arc_Q*obj.arc.sol_s(stp,1);

    end

    function comp_stp(obj)
    
      if obj.dgnflg
        obj.stp = obj.stpmax;
        obj.x = obj.x + obj.stp*obj.p;
        obj.info_srch = 0;
        if obj.opt.print_dbg
          obj.dprint('search','stp','%9.2e',obj.stp);
        end
        return;
      end
      
      switch obj.opt.search
        case 'line'
          arc_handle = @(stp_param) line_arc(stp_param,obj.p);
        case 'arc'
          arc_handle = @(stp_param) obj.arc_eval(stp_param);
      end
      
      fevmax = min(obj.opt.fevmax_srch,obj.opt.fevmax-obj.fevcnt);
      ncur = 0;
      [obj.x obj.f obj.g obj.stp obj.info_srch fev_srch] = mtasrch(...
        obj.usrfcn,obj.x,obj.f,obj.g,arc_handle,obj.stpinit,...
        ncur,obj.opt.ftol,obj.opt.gtol,obj.opt.xtol,...
        obj.opt.stpmin,obj.stpmax,fevmax,0,obj.opt.bisect_srch);
      
      if obj.stpbnd == 0 || ...
          ( obj.info_srch ~= 2 && obj.info_srch ~= 5 && obj.info_srch ~= 6 )
        % a bound was not hit
        obj.stpbnd = 0;
        obj.stpidx = 0;
      end

      % update function eval counter
      obj.fevcnt = obj.fevcnt + fev_srch;
      
      if obj.opt.print_dbg
        obj.dprint('search','stp','%9.2e',obj.stp);
        obj.dprint('search','info_srch','%d',obj.info_srch);
        obj.dprint('search','fev_srch','%d',fev_srch);
        
%         stp = obj.stp;
%         fp = obj.f;
%         finit = obj.f_old;
%         ginit = obj.p'*obj.g_old;
%         dp = obj.p'*obj.g;
%         if fp <= finit + obj.opt.ftol*(ginit*stp + 0.5*min(ncur,0)*stp^2) ...
%             || abs(dp) <= obj.opt.gtol*abs(ginit + min(ncur,0)*stp)
%           obj.dprint('search','srch cond satisfied','%s','yes');
%         else
%           obj.dprint('search','srch cond satisfied','%s','no');
%         end
        
      end
      
    end

    function update_B(obj)
      
      update_flg = 0;
      if strcmp(obj.opt.update,'BFGS') && ismember(obj.info_srch,[1 6]) && ~obj.dgnflg
        s_up = obj.x(obj.xfree)-obj.x_old(obj.xfree);
        y_up = obj.g(obj.xfree)-obj.g_old(obj.xfree);
        Bs = obj.B(obj.xfree,obj.xfree)*s_up;
        obj.B(obj.xfree,obj.xfree) = obj.B(obj.xfree,obj.xfree) - Bs*Bs'/(s_up'*Bs) + y_up*y_up'/(y_up'*s_up);
        update_flg = 1;
      elseif strcmp(obj.opt.update,'SR1') && ismember(obj.info_srch,[1 5 6]) && ~obj.dgnflg
        s_up = obj.x-obj.x_old;
        y_up = obj.g-obj.g_old;
        ymBs = y_up-obj.B*s_up;
        
        if abs(s_up'*ymBs) > obj.opt.sr1_tol*norm(s_up)*norm(ymBs)
          obj.B = obj.B + (ymBs*ymBs') / (s_up'*ymBs);
          update_flg = 1;
        end
        
      end
      
      % clear out old bfgs info when a constraint is hit
      delete_flg = 0;
      if obj.stpidx && strcmp(obj.opt.search,'line')
        obj.B(obj.stpidx,:) = 0;
        obj.B(:,obj.stpidx) = 0;
        obj.B(obj.stpidx,obj.stpidx) = 1;
        delete_flg = 1;
      end
      
      if obj.opt.print_dbg
        cond_B = cond(obj.B);
        evals = eig(obj.B);
        symm_check = max(max(abs(obj.B-obj.B')));
        
        obj.dprint('update_B','update_flag','%d',update_flg);
        obj.dprint('update_B','delete_flag','%d',delete_flg);
        obj.dprint('update_B','cond(B)','%8.2e',cond_B);
        obj.dprint('update_B','min(eig(B))','%8.2e',min(evals));
        obj.dprint('update_B','max(eig(B))','%8.2e',max(evals));
        obj.dprint('update_B','norm(B-Bt)','%8.2e',symm_check);
        
        
        if min(evals) <= 0 || cond_B >= obj.opt.infval || symm_check >= obj.opt.ptol
          obj.dprint('update_B','exception','%s','bad B');
        end
        
      end
      
    end
    
    function update_state(obj)
      
      % if a non-degenerate step is taken, unfix variables
      if ~obj.dgnflg || ~any(obj.xstate == 0)
        obj.xstate(obj.xstate == -2) = -1;
        obj.xstate(obj.xstate == 2) = 1;
      end
      
      % deactivate bound hit in this iteration
      if obj.stpbnd == -1
        obj.xstate(obj.stpidx) = -2;
      elseif obj.stpbnd == 1
        obj.xstate(obj.stpidx) = 2;
      end
        
      % free variables if appropriate
      obj.xstate(obj.xstate == -1 & obj.g < 0) = 0;
      obj.xstate(obj.xstate == 1 & obj.g > 0) = 0;
      
      % set xfree
      obj.xfree = obj.xstate == 0;
      
      if obj.opt.print_dbg
        obj.dprint('update_state','sum(xstate==-2)','%d',sum(obj.xstate==-2));
        obj.dprint('update_state','sum(xstate==-1)','%d',sum(obj.xstate==-1));
        obj.dprint('update_state','sum(xstate==1)','%d',sum(obj.xstate==1));
        obj.dprint('update_state','sum(xstate==2)','%d',sum(obj.xstate==2));
        obj.dprint('update_state','sum(xfree)','%d',sum(obj.xfree));
      end
      
    end
    
    %% expand methods
    function exp_init(obj)
      %expand_init  initialize expand values
      %
      % On entry:
      %   opt.ptol = primal tolerance
      %   opt.expa = initial feasibility parameter
      %   opt.expb = final feasibility parameter
      %   opt.expfrq = expand reset frequency
      %
      % On exit:
      %   exptol = initial expand tolerance
      %   expinc = expand increase value
      %   expcnt = 0, because expand has just been initialized
      %
      
      obj.exptol = obj.opt.expa*obj.opt.ptol;
      obj.expinc = (obj.opt.expb - obj.opt.expa)*obj.opt.ptol / obj.opt.expfrq;
      obj.expcnt = 0;
      
    end
    
    function expflg = exp_main(obj)
      %exp_main  expand or reset working feasibility tolerance
      %
      % Method calls:
      %   exp_reset
      %
      % On entry:
      %   expcnt = number of expand increases since last reset
      %   opt.expfrq = number of allowed expand increases before reset
      %   exptol = current expand feasibility tolerance
      %   expinc = expand increase value
      %
      % On exit:
      %   expcnt = incremented by 1 or set to 0 if reset
      %   exptol = increased or reset feasibility tolerance
      %
      % Returns:
      %   expflg = true if expand reset results in infeasibility
      %            false otherwise
      %
      
      expflg = false;
      
      if obj.expcnt < obj.opt.expfrq
        % expand the working feasibility tolerance
        obj.exptol = obj.exptol + obj.expinc;
        obj.expcnt = obj.expcnt + 1;
      else
        % must reset
        expflg = obj.exp_reset();
      end
      
    end
    
    function expflg = exp_reset(obj)
      %exp_reset  reset procedure for expand
      %
      % Returns:
      %   expflg = true if reset results in an infeasible point
      %            false if reset produces a feasible point
      %
      
      % TODO fix this
      
      expflg = false;
      
      % reinitialize expand
      obj.exp_init();
      
    end
    
    
  end
  
end









































