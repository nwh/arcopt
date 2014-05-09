classdef arcopt_nm_lc < handle
  %arcopt_nm_lc ARCOPT: Newton Method for problems with Linear Constraints
  %

  properties

    %% general options
    opt = []; % options structure

    %% problem size
    n = []; % number of variables
    m = []; % number of constraints

    %% workspace
    wn = []; % workspace vector of size n
    wm = []; % workspace vector of size m
    wmn = []; % workspace vector of size m+n

    %% input data
    usrfun = [];  % user function handle
    usrhess = []; % user hessian handle
    x0 = [];     % initial point
    A_in = [];   % constraint matrix
    bl_in = [];  % lower bound on variable and constraint
    bu_in = [];  % upper bound on variable and constraint
    noeq_mode = []; % flag is true if there are no linear equality constraints

    %% working variables
    x = [];  % current iterate, length = n+m, contains slacks
    x_old =[]; % previous iterate
    f = [];  % working function value
    f_old = []; % previous function value
    g = [];  % working gradient vector
    rg = []; % working reduced gradient vector
    H = [];  % working Hessian matrix
    p = []; % (modified) newton direction
    bl = []; % working lower bound
    bu = []; % working upper bound
    A = [];  % full matrix, [A_in -I]
    y = [];  % multiplier vector
    z = [];  % residual gradient vector
    bl_lx = []; % logical index indicating bound violation
    bu_lx = []; % logical index indicating bound violation

    dx = []; % the search direction
    rdx = []; % the reduced search direction
    sigma = []; % sign of free variable in dx

    %% dnc, negative curvature
    dnc = []; % direction of negative curvature in reduced space
    dncflg = 0; % negative curvature flag.  true if dnc exists, false otherwise
    dnccnt = 0; % dnc counter
    dncmflg = []; % dnc modification flag
                  % true if dnc is added to arc rhs
                  % false otherwise
    dncmcnt = 0; % arc modification counter

    %% basis factorization
    bfac_opt = []; % basis factorization options
    bfac = []; % basis factorization object
    bcol = []; % basis column list
               % A(:,bcol(i)) is in column i of B
    depcol = []; % logical array indicating dependent columns in B
    nsing = []; % number of apparent singularities in B
    s_ix = []; % index of variable to become superbasic
    r_ix = []; % index of variable to become nonbasic
    b_ix = []; % index of variable to become basic
    hs = []; % basis information vector
             % hs(i) = -1 if x(i) is nonbasic between bl(i) and bu(i)
             % hs(i) = 0 if x(i) is nonbasic at bl(i)
             % hs(i) = 1 if x(i) is nonbasic at bu(i)
             % hs(i) = 2 if x(i) is superbasic
             % hs(i) = 3 if x(i) is basic
             % hs(i) = 4 if x(i) is nonbasic and fixed (bl(i) == bu(i))
    facucnt = 0; % factorization update counter

    %% factorization flags
    upflg = false; % true if basis update attempted
    reflg = false; % true if basis refactorization attempted
    bsflg = false; % true if BS factorization called
    brflg = false; % true if BR factorization called

    %% termination
    termchk = []; % value for termination check

    %% step size properties
    stp = []; % the search step
    stpinit = []; % initial step size
    stpmax = []; % step to a bound or to opt.stpmax
    stpbnd = []; % -1 if bl hit, 1 if bu hit, 0 if no bound hit

    %% arc properties
    arc = []; % arc object
    arc_dim = []; % arc dimension
    arc_Q = []; % arc subspace matrix
    arc_ZQ = []; % arc nullspace matrix
    arc_H = []; % arc matrix term
    arc_g = []; % arc inhomogeneous term

    %% eigs
    opt_eigs = []; % options structure for eigs
    eval = []; % minimum eigenvalue of reduced Hessian
    evalflg = []; % true if negative eigenvalue is present
    evec = []; % eigenvector of reduced Hessian corresponding to eval
    eigsflg = []; % output flag for eigs

    %% linear system solve properties
    linflg = []; % lin flag
    lincnt = []; % lin iteration count
    linres = []; % lin relative residual

    %% expand variables
    expcnt = 0; % expand counter
    exptol = 0; % working feasibility tolerance
    expinc = 0; % expand increase value (tau in paper)

    %% degeneracy
    dgnflg = 0; % degenerate step flag
    dgn1cnt = 0; % phase 1 degenerate step counter
    dgn2cnt = 0; % phase 2 degenerate step counter

    %% arc search properties
    srch_flg = 0; % status flag from arc search
    srch_fev = 0; % number of function evaluations from arc search
    srch_d2val = 0; % initial value for second derivative of search function

    %% counters
    fevcnt = 0; % function evaluation counter
    prntcnt = 1; % print counter
    itercnt = 0; % iteration counter

    %% phase 1 data
    infcnt = 0; % infeasibility count
    infsum = 0; % sum of infeasibilities
    phase1cnt = 0; % phase 1 iteration counter

    %% phase 2 data
    phase2cnt = 0; % phase 2 iteration counter

  end

  methods (Static)

    function options = optset(varargin)
      %optset  construct options structure with default vaules
      %
      %

      % get the input parser
      in_parse = inputParser;

      % basic parameters
      in_parse.addParamValue('dtol',1e-6,@(x) x>=0); % dual convergence tolerance
      in_parse.addParamValue('ptol',1e-6,@(x) x>=0); % primal convergence tolerance
      in_parse.addParamValue('ctol',1e-4,@(x) x>=0); % curvature convergence tolerance
      in_parse.addParamValue('itermax',1000,@(x) x>0); % maximum number of iterations
      in_parse.addParamValue('fevmax',1000,@(x) x>0); % maxumum number of function evaluations
      in_parse.addParamValue('facfrq',50,@(x) x>=0); % refactorization frequency
      in_parse.addParamValue('infval',1e20,@(x) x>0); % infinity value
      in_parse.addParamValue('cycle_guard',0,@(x) x==0 || x==1 || islogical(x)); % cycle guard
      
      % crash procedure
      in_parse.addParamValue('crash','slack',@(x) ismember(x,{'slack','firstm'})); % selector for crash prodecure
      
      % expand parameters
      in_parse.addParamValue('expfrq',10000,@(x) x>=0); % expand reset frequency
      in_parse.addParamValue('expa',.5,@(x) x>0 && x<1); % parameter for initial feasibility tolerance
      in_parse.addParamValue('expb',.99,@(x) x>0 && x<1); % parameter for final feasibility tolerance
      in_parse.addParamValue('expsml',1e-11,@(x) x>=0); % tolerance for near zero numbers

      % arc search parameters
      in_parse.addParamValue('ftol',1e-4,@(x) x>0 && x<1); % arc search descent parameter
      in_parse.addParamValue('gtol',.9,@(x) x>0 && x<1); % arc search curvature parameter
      in_parse.addParamValue('xtol',1e-8,@(x) x>0); % arc search interval tolerance
      in_parse.addParamValue('fevas',100,@(x) x>0); % arc search max function evaluations
      in_parse.addParamValue('stpmax',1e10,@(x) x>0); % arc search max step size
      in_parse.addParamValue('stpmin',0,@(x) x>=0); % arc search min step size

      % dnc options
      in_parse.addParamValue('dnc_mtol',0.2,@(x) 0<x && x<1.3); % arc modification tolerance

      % arc parameters
      in_parse.addParamValue('arc_vtol',1e-7,@(x) x>=0); % tolerance for intial step size selection

      % eigs options
      in_parse.addParamValue('eigs_tol',sqrt(eps),@(x) x>0); % tolerance for eigs convergence
      in_parse.addParamValue('eigs_maxit',300,@(x) x>0); % maximum number of iterations for eigs
      in_parse.addParamValue('eigs_p',20,@(x) x>0); % number of vectors used in eigs
      in_parse.addParamValue('eigs_disp',0,@(x) ismember(x,[0 1 2])); % eigs display option

      % linear system solve options
      in_parse.addParamValue('lin_tol',1e-6,@(x) x>0); % tolerance for linear system solve convergence
      in_parse.addParamValue('lin_maxit',[]); % maximum number of iterations for each linear system solve
      in_parse.addParamValue('lin_delta',1e-4,@(x) x>=0); % shifting factor for linear solves

      % print options
      in_parse.addParamValue('print_level','iter',@(x) ismember(x,{'iter','stats','none'}));
      in_parse.addParamValue('print_screen',1,@(x) ismember(x,[0 1]));
      in_parse.addParamValue('print_file',0,@(x) x == 0 | x > 2);
      in_parse.addParamValue('print_freq',10,@(x) x>0);

      % parse the input
      in_parse.parse(varargin{:});

      % obtain the output
      options = in_parse.Results;

    end

  end

  methods

    %% constructor & options setup
    function obj = arcopt_nm_lc(usrfun,usrhess,x0,bl,bu,A,cl,cu,varargin)

      % process options
      obj.opt = obj.optset(varargin{:});

      % process eigs options
      obj.eigs_options();

      % get number of variables from initial point x0
      if isvector(x0)
        obj.n = length(x0);
      else
        error('arcopt_nm_lc:input','x0 must be a vector.');
      end

      % if A is not specified
      if isempty(A)
        obj.mprint('stats','\nAdding dummy linear equality.\n')
        obj.noeq_mode = true;
        obj.m = 0;
        n_temp = obj.n;
      else
        obj.noeq_mode = false;
        [obj.m n_temp] = size(A);
      end

      % check input sizes

      % check size of A
      if n_temp ~= obj.n
        error('arcopt_nm_lc:input','A must have n columns.')
      end

      if length(bl) ~= obj.n
        error('arcopt_nm_lc:input','bl must have length n.');
      end

      if length(bu) ~= obj.n
        error('arcopt_nm_lc:input','bu must have length n.');
      end

      if length(cl) ~= obj.m
        error('arcopt_nm_lc:input','cl must have length m.');
      end

      if length(cu) ~= obj.m
        error('arcopt_nm_lc:input','cu must have length m.');
      end

      % modify input data if there are no equality constraints
      if obj.noeq_mode
        A = zeros(1,obj.n);
        cl = -obj.opt.infval;
        cu = obj.opt.infval;
        obj.m = 1;
      end

      % store in object
      obj.usrfun = usrfun;
      obj.usrhess = usrhess;
      obj.x0 = x0(:);
      obj.A_in = A;
      obj.bl_in = [bl(:); cl(:)];
      obj.bu_in = [bu(:); cu(:)];

    end

    function eigs_options(obj)
      %eigs_options  set options for eigs
      
      obj.opt_eigs.issym = 1;
      obj.opt_eigs.isreal = 1;
      obj.opt_eigs.tol = obj.opt.eigs_tol;
      obj.opt_eigs.maxit = obj.opt.eigs_maxit;
      obj.opt_eigs.p = obj.opt.eigs_p;
      obj.opt_eigs.disp = obj.opt.eigs_disp;
      
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
        error('arcopt_nm_lc:mprint','invalid print level string: %s.',level);
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

    function print_start(obj)
      %print_start print initial stats

      obj.mprint('stats','\n')

      obj.mprint('stats','Matrix stats:\n');
      obj.mprint('stats',' m = %d\n',obj.m);
      obj.mprint('stats',' n = %d\n',obj.n);
      obj.mprint('stats',' nnz(A) = %d\n',nnz(obj.A_in))
      obj.mprint('stats',' max(|A|) = %7.1e\n',full(max(max(abs(obj.A_in)))));

      obj.mprint('stats','\n')

      % compute bound status for printing
      % bs(i) = 0 if bl(i) > -inf and bu(i) = inf
      % bs(i) = 1 if bl(i) = -inf and bu(i) < inf
      % bs(i) = 2 if bl(i) > -inf and bu(i) < inf and |bl(i)-bu(i)| > ptol
      % bs(i) = 3 if bl(i) > -inf and bu(i) < inf and |bl(i)-bu(i)| <= ptol
      % bs(i) = 4 if bl(i) = -inf and bu(i) = inf
      bs = zeros(obj.n+obj.m,1);
      for i = 1:(obj.n+obj.m)
        if obj.bl(i) > -obj.opt.infval && obj.bu(i) >= obj.opt.infval
          % x(i) is bdd below
          bs(i) = 0;
        elseif obj.bl(i) <= -obj.opt.infval && obj.bu(i) < obj.opt.infval
          % x(i) is bdd above
          bs(i) = 1;
        elseif abs(obj.bl(i)-obj.bu(i)) > obj.opt.ptol
          % x(i) is bdd above and below
          bs(i) = 2;
        elseif abs(obj.bl(i)-obj.bu(i)) <= obj.opt.ptol
          % x(i) is fixed
          bs(i) = 3;
        elseif obj.bl(i) <= -obj.opt.infval && obj.bu(i) >= obj.opt.infval
          % x(i) is unbounded
          bs(i) = 4;
        else
          % this should not happen!
          error('arcopt_nm_lc:start','bad bounds.')
        end
      end

      obj.mprint('stats','Linear constraints:\n')
      obj.mprint('stats',' <- : %d\n',sum(bs(obj.n+1:end)==0));
      obj.mprint('stats',' -< : %d\n',sum(bs(obj.n+1:end)==1));
      obj.mprint('stats',' << : %d\n',sum(bs(obj.n+1:end)==2));
      obj.mprint('stats',' == : %d\n',sum(bs(obj.n+1:end)==3));
      obj.mprint('stats',' -- : %d\n',sum(bs(obj.n+1:end)==4));

      obj.mprint('stats','\n')

      obj.mprint('stats','Variable bounds:\n')
      obj.mprint('stats',' <- : %d\n',sum(bs(1:obj.n)==0));
      obj.mprint('stats',' -< : %d\n',sum(bs(1:obj.n)==1));
      obj.mprint('stats',' << : %d\n',sum(bs(1:obj.n)==2));
      obj.mprint('stats',' == : %d\n',sum(bs(1:obj.n)==3));
      obj.mprint('stats',' -- : %d\n',sum(bs(1:obj.n)==4));

      %obj.mprint('stats','\n')

    end

    %% main solve call and initialization methods
    function [x f flg aux] = solve(obj)
      %solve  solve the problem

      % initiate timer
      solvetime = tic;

      % initialize counters
      obj.phase1cnt = 0;
      obj.phase2cnt = 0;
      obj.itercnt = 0;

      % initialize solve flag
      % everything is good at this point
      solveflg = 1;

      % augment the matrix
      % the last m columns correspond to the slack variables
      obj.A = [obj.A_in -eye(obj.m)];

      % initialize workspace
      obj.init_workspace();

      % check feasibility of bounds
      infbnd = obj.bnd_feas_chk();
      % obj.bl and obj.bu are now set

      if infbnd > 0
        % bounds are infeasible
        solveflg = 11;
      end

      % initial print
      obj.print_start();

      if solveflg == 1
        % load initial point
        obj.init_x();

        % construct intial hs
        switch obj.opt.crash
          case 'slack'
            obj.basis_init_slack();
          case 'firstm'
            obj.basis_init_firstm();
        end

        % initialize expand
        obj.exp_init();

        % initialize basis factorization
        obj.fac_init();
        obj.fac_main();

        % initiate phase1
        solveflg = obj.phase1();

      end

      if solveflg == 101
        % a feasible point has been found, begin phase 2
        solveflg = obj.phase2();
      end

      % stop timer
      solvetime = toc(solvetime);

      % set final variables
      x = obj.x(1:obj.n);
      f = obj.f;
      flg = solveflg;

      if obj.noeq_mode
        aux.s = [];
        aux.y = [];
      else
        aux.s = obj.x(obj.n+1:end);
        aux.y = obj.y;
      end

      % set up auxiliary output structure
      aux.g = obj.g(1:obj.n);
      aux.z = obj.z;
      aux.hs = obj.hs;
      aux.phase1cnt = obj.phase1cnt;
      aux.phase2cnt = obj.phase2cnt;
      aux.dng1cnt = obj.dgn1cnt;
      aux.dng2cnt = obj.dgn2cnt;
      aux.dnccnt = obj.dnccnt;
      aux.dncmcnt = obj.dncmcnt;
      aux.itercnt = obj.itercnt;
      aux.fevcnt = obj.fevcnt;
      aux.time = solvetime;
      aux.srch_flg = obj.srch_flg;

      % print final stats
      obj.mprint('stats','\n');
      obj.mprint('stats','Termination status: ');
      switch solveflg

        % exit codes from solve method
        case 11
          obj.mprint('stats','infeasible bounds.\n')

          % exit codes from phase 1
        case 111
          obj.mprint('stats','phase 1 reports infeasible constraints.\n')
        case 131
          obj.mprint('stats','phase 1 exhausted maximum iterations.\n')

          % exit codes from phase 2
        case 1
          obj.mprint('stats','optimality conditions satisfied.\n')
        case 3
          obj.mprint('stats','phase 2 exhausted maximum iterations.\n')
        case 2
          obj.mprint('stats','phase 2 exhausted maximum function evaluations.\n')

        otherwise
          obj.mprint('stats','unknown.\n')
      end

      obj.mprint('stats','Function value....... %g\n',obj.f);
      obj.mprint('stats','Residual............. %g\n',norm(obj.A*obj.x,'inf'))
      obj.mprint('stats','Phase 1 iterations... %d\n',obj.phase1cnt);
      obj.mprint('stats','Phase 2 iterations... %d\n',obj.phase2cnt);
      obj.mprint('stats','Total iterations..... %d\n',obj.itercnt);
      obj.mprint('stats','dnc iterations....... %d\n',obj.dnccnt);
      obj.mprint('stats','dnc modifications.... %d\n',obj.dncmcnt);
      obj.mprint('stats','P1 degenerate steps.. %d\n',obj.dgn1cnt);
      obj.mprint('stats','P2 degenerate steps.. %d\n',obj.dgn2cnt);
      obj.mprint('stats','Function evaluations. %d\n',obj.fevcnt);
      obj.mprint('stats','Time (sec)........... %g\n',solvetime);

    end

    function init_workspace(obj)
      %init_workspace  initialize workspace vectors to correct size

      obj.bcol = zeros(obj.m,1);
      obj.depcol = false(obj.m,1);

      obj.wm = zeros(obj.m,1);
      obj.wn = zeros(obj.n,1);
      obj.wmn = zeros(obj.m+obj.m,1);

      obj.f = 0;
      obj.x = zeros(obj.n+obj.m,1);
      obj.hs = zeros(obj.n+obj.m,1);
      obj.g = zeros(obj.n+obj.m,1);
      obj.bl = zeros(obj.n+obj.m,1);
      obj.bu = zeros(obj.n+obj.m,1);
      obj.z = zeros(obj.n+obj.m,1);
      obj.y = zeros(obj.m,1);

      obj.bl_lx = false(obj.n+obj.m,1);
      obj.bu_lx = false(obj.n+obj.m,1);

      obj.dx = zeros(obj.n+obj.m,1);
      obj.stp = 0;

    end

    function infbnd = bnd_feas_chk(obj)
      %bnd_feas_chk  simple feasibility check on bounds
      %
      % On entry:
      %   bl_in = input lower bounds
      %   bu_in = input upper bounds
      %   opt.infval = infinity value
      %
      % On exit:
      %   bl_in = min val is now -opt.infval
      %   bu_in = max val is now opt.infval
      %   bl = working lower bound
      %   bu = working upper bound
      %
      % Output:
      %   infbnd = 1 if bounds are infeasible
      %            0 if bounds are feasible
      %

      % make sure bounds are not outside [-infval infval]
      obj.bl_in(obj.bl_in < -obj.opt.infval) = -obj.opt.infval;
      obj.bu_in(obj.bu_in > obj.opt.infval) = obj.opt.infval;

      % set working bounds
      obj.bl = obj.bl_in;
      obj.bu = obj.bu_in;

      infbnd = 0;

      % no lower bound can be >= infval
      if any(obj.bl >= obj.opt.infval)
        infbnd = 1;
      end

      % no upper bound can be <= -infval
      if any(obj.bu <= -obj.opt.infval)
        infbnd = 1;
      end

      % no lower bound can be greater than an upper bound
      if any(obj.bl > obj.bu)
        infbnd = 1;
      end

    end

    function init_x(obj)
      %init_x  intialize point and basis information
      %
      % Initialized the variable vector x.  It takes the user input x0 and
      % projects it into the bounds.
      %
      % On entry:
      %   x = array of size n+m
      %   x0 = input intial point
      %   n = number of variables
      %   bl = initialized lower bound
      %   bu = initialized upper bound
      %
      % On exit:
      %   x(1:n) = initial point projected into bounds
      %   x(n+1:m) = slack variables, all zero
      %

      % load initial point
      obj.x(:) = 0;
      obj.x(1:obj.n) = obj.x0;
      obj.x(obj.n+1:obj.n+obj.m) = obj.A_in*obj.x0;
      if any(obj.x < obj.bl) || any(obj.x > obj.bu)
        % project x into bounds
        obj.mprint('stats','Initial point projected into bounds.\n');
        obj.x = max(obj.x,obj.bl);
        obj.x = min(obj.x,obj.bu);
      end

    end

    %% wrappers for user function calls
    function [f g] = func_wrap(obj,x)
      %func_wrap  wrapper for user function call
      g = zeros(obj.n+obj.m,1);
      [f g(1:obj.n)] = obj.usrfun(x(1:obj.n));
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
      % Method calls:
      %   restore
      %   fac_main
      %   exp_init
      %   comp_basics
      %   basic_feas_chk
      %
      % On entry:
      %   !! many methods called
      %
      % On exit:
      %   !! many methods called
      %
      % Returns:
      %   expflg = true if reset results in an infeasible point
      %            false if reset produces a feasible point
      %

      expflg = false;

      % restore infeasible superbasics and nonbasics to their bounds
      obj.restore();

      % refactorize the matrix
      obj.fac_main();

      % reinitialize expand
      obj.exp_init();

      % recompute basic variables
      obj.comp_basics();

      % check feasibility of basics
      obj.basic_feas_chk();

      % if resulting point is infeasible, return expflg as true
      if obj.infcnt > 0
        expflg = true;
      end

    end

    %% factorization methods
    function fac_init(obj)
      %fac_init  instantiate the basis factorization object
      %
      % On exit:
      %   bfac_opt = basis factorization options structure
      %   bfac = empty basis factorization object
      %

      % (TODO) expose lusol options to user
      obj.bfac_opt = lusol.luset();
      obj.bfac = lusol(1,obj.bfac_opt);

    end

    function fac_main(obj)
      %fac_main  handle factorization and basis repair
      %
      % This method carries out LU factorization of the basis matrix.  If
      % the factorization indicates singularity, then basis repair is
      % attempted.
      %
      % The general flow of the procedure:
      %   1) factorize B
      %   2) if singular, envoke BS factorization.  This will swap columns
      %      of B and S to achieve non-singular and better conditioned B.
      %   3) if still singular, envoke BR factorization.  This will replace
      %      dependent columns of B with appropriate columns of I
      %      corresponding to slack variables.  This will certainly produce
      %      a non-singular basis, because all columns could be replaced.
      %   4) If (2) and/or (3) occured, refactorize B for use in algorithm.
      %
      % Possible method calls:
      %   fac_BS
      %   fac_BR
      %   fac_B
      %
      % On entry:
      %   A = constraint matrix with explicit columns for slack variables
      %   hs = basis information vector
      %   bfac = basis factorization object
      %   bcol = basis information vector
      %
      % On exit:
      %   hs = may have changed because of BS or BR factorization
      %   bfac = contains non-singular factorized matrix
      %   bcol = basis column map
      %   reflg = true, because (re)factorization happened
      %   bsflg = true if BS factorization called
      %   brflg = true if BR factorization called
      %

      % set refactorization flag
      obj.reflg = true;
      obj.bsflg = false;
      obj.brflg = false;

      % attempt to factorize B
      obj.fac_B();

      if obj.nsing > 0 && sum(obj.hs == 2) > 0
        % basis appears singular, there are superbasic variables
        % attempt BS factorization
        obj.fac_BS();
        obj.bsflg = true;
      end

      if obj.nsing > 0
        % basis appears singular, there are no superbasic variables
        % attempt BR factorization
        obj.fac_BR();
        obj.brflg = true;
      end

      if obj.bsflg || obj.brflg
        % Need to perform a full factorization after BS or BR
        % factorization.
        obj.fac_B();
      end

      % reset factorization update counter
      obj.facucnt = 0;

    end

    function fac_B(obj)
      %fac_B  factorize basis matrix
      %
      % Here we simply factorize B.  LUSOL is instructed to use threshold
      % partial pivoting and keep the LU factors for use in solves.
      %
      % Method calls:
      %   bfac.factorize
      %
      % On entry:
      %   hs = basis information vector
      %        here we assume that sum(hs == 3) == m
      %
      % On exit:
      %   bfac = basis factorization object ready for solves
      %   bcol = basis column map
      %   depcol = indication of dependent columns
      %            this will likely not be used for basis repair
      %            the basis repair method will use the TRP option
      %   nsing = number of aparent singularities
      %

      obj.bfac_opt.pivot = 'TPP';
      obj.bfac_opt.keepLU = 1;
      obj.bcol = find(obj.hs == 3);
      [inform obj.nsing obj.depcol] = ...
        obj.bfac.factorize(obj.A(:,obj.bcol),obj.bfac_opt);

    end

    function fac_BS(obj)
      %fac_BR  swap columns of B and S to obtain for better cond(B)
      %
      % Say B = A(:,hs==3) is the matrix of basic columns
      % Say S = A(:,hs==2) is the matrix of superbasic columns
      %
      % This method performs a factorization on [B S]' to find a new set of
      % columns that will make a better conditioned B.
      %
      % LUSOL is told to use threshold rook pivoting, which helps to
      % identify dependent columns.  The LU factors are thrown away,
      % because they will not be used for solves.
      %
      % Method calls:
      %   bfac.factorize
      %
      % On entry:
      %   hs = basis information vector
      %        we assume there are m basic variables
      %        we assume there are more than 1 superbasic variables
      %
      % On exit:
      %   hs = updated basis information vector
      %   nsing = number of apparent singularities
      %   bfac = has been modified, cannot be used for solves!
      %   bfac_opt = has been modified
      %

      % set lusol options
      obj.bfac_opt.pivot = 'TRP';
      obj.bfac_opt.keepLU = 0;

      % get indecies basic and superbasic columns
      cols = find(obj.hs == 3 | obj.hs == 2);

      % number of columns in this factorization
      nbs = length(cols);

      % perform the factorization
      [inform obj.nsing obj.depcol] = ...
        obj.bfac.factorize(obj.A(:,cols)',obj.bfac_opt);

      % update the basis information vector
      obj.hs(cols(obj.bfac.ip(1:obj.m))) = 3;
      obj.hs(cols(obj.bfac.ip(obj.m+1:nbs))) = 2;

      % clear out basis column map
      obj.bcol(:) = 0;

    end

    function fac_BR(obj)
      %fac_BR  replace columns of B with slacks to obtain full rank basis
      %
      % The BR factorization uses LUSOL's threshold rook pivoting option to
      % find dependent columns.  Dependent columns are replaced with
      % identity columns corresponding to slack variables in the fac_repair
      % method.
      %
      % There is no guarantee that one factorization and subsequent call to
      % fac_repair will produce a nonsingular basis.  Thus, this method
      % will repeat until a full rank basis is produced.  This is always
      % possible.  Just imagine the case where the entire basis matrix was
      % replaced with the identity.
      %
      % For efficiency the method throws the LU factors away.  A
      % subsequent factorization is needed for solves.
      %
      % Method calls:
      %   bfac.factorize
      %   fac_repair
      %
      % On entry:
      %   hs = basis information vector
      %   bfac = basis factorization object
      %   nsing = apparent number of singularities.  Basis repair happens
      %           if nsing >= 1.
      %
      % On exit:
      %   bfac = modified basis factorization class.  Not ready for solves.
      %   bfac_opt = modified basis factorization options
      %   hs = modified basis information vector giving full rank basis
      %

      % set lusol options
      obj.bfac_opt.pivot = 'TRP';
      obj.bfac_opt.keepLU = 0;

      % force basis repair check
      obj.nsing = 1;

      while( obj.nsing )
        % perform BR factorization to find dependent columns
        obj.bcol = find(obj.hs == 3);
        [inform obj.nsing obj.depcol] = obj.bfac.factorize( ...
          obj.A(:,obj.bcol),obj.bfac_opt);

        if obj.nsing
          % dependent columns found, replace with slack variables
          obj.fac_repair();
        end
      end

    end

    function fac_repair(obj)
      %fac_repair  factorization repair, adapted from s2sing in SNOPT
      %
      % This method carries out factorization or basis repair by making:
      %   variables corresponding to dependent columns nonbasic
      %   appropriate slack variables basic
      %
      % Before factorizing matrix A, LUSOL searches for all rows and
      % columns of length 1 and permutes them to the top-left.  If there
      % are duplicate singleton rows or columns, only one of them is
      % permuted.  When this is complete, elimination begins on the
      % bottom-right.
      %
      % Basis repair operates by replacing dependent columns of B with
      % singleton columns corresponding to slack variables.
      %
      % LUSOL returns information on which columns it determines are
      % dependent, this is stored in logical vector depcol.  This method is
      % based on SNOPT's s2sing and swaps variables to result in a full
      % rank basis.
      %
      % On entry:
      %   bfac = basis factorization object, this method uses the
      %          permutation vectors ip and iq
      %   bcol = basis column map
      %   hs = basis information vector
      %   x = variable vector
      %
      % On exit:
      %   hs = basis variables are swapped with appropriate slacks
      %   x = some basic variables may have been moved and placed on their
      %       bounds
      %

      for k = 1:obj.m
        j = obj.bfac.iq(k);
        if obj.depcol(j)
          j = obj.bcol(j);
          % column corresponding to variable j is dependent
          % make j nonbasic
          if obj.x(j) <= obj.bl(j)
            obj.x(j) = obj.bl(j);
            obj.hs(j) = 0;
          elseif obj.x(j) >= obj.bu(j)
            obj.x(j) = obj.bu(j);
            obj.hs(j) = 1;
          else
            obj.hs(j) = -1;
          end

          if obj.bl(j) == obj.bu(j)
            % variable j is fixed
            obj.hs(j) = 4;
          end

          i = obj.bfac.ip(k);
          % make i basic
          obj.hs(obj.n+i) = 3;
        end
      end
    end

    %% basis methods
    function basis_init_slack(obj)
      %basis_init_slack  basis is initialized with all slack variables
      %
      % Initialize the basis.  This simply selects all slack variables as being
      % in the basis.  This should be replaced by a better crash procedure.
      % This method does make the first factorization trivial to factorize.
      %
      % Any variable that is initially at an upper or lower bound is set as
      % nonbasic.  This may not be the best thing to do here.
      %
      % On entry:
      %   x = initialized variable vector
      %   bl = initialized lower bound
      %   bu = initialized upper bound
      %   n = number of variables
      %   hs = basis information vector, size n+m
      %
      % On exit:
      %   hs = initialized basis vector, all slacks are basic
      %

      % everything is superbasic initially
      obj.hs(:) = 2;

      % make variables at lower bound nonbasic
      obj.hs( obj.x == obj.bl ) = 0;

      % make variables at upper bound non basic
      obj.hs( obj.x == obj.bu ) = 1;

      % set fixed variables
      obj.hs( obj.bl == obj.bu ) = 4;

      % all slacks are basic
      obj.hs(obj.n+1:end) = 3;

    end

    function basis_init_firstm(obj)
      %basis_init_firstm  basis is initialized with first m variables
      %
      % Initialize the basis.  This simply selects the first m variables as being
      % in the basis.  This should be replaced by a better crash procedure.
      %
      % Any variable that is initially at an upper or lower bound is set as
      % nonbasic.  This may not be the best thing to do here.
      %
      % On entry:
      %   x = initialized variable vector
      %   bl = initialized lower bound
      %   bu = initialized upper bound
      %   m = number of constraints
      %   hs = basis information vector, size n+m
      %
      % On exit:
      %   hs = initialized basis vector, all slacks are basic
      %

      % everything is superbasic initially
      obj.hs(:) = 2;

      % make variables at lower bound nonbasic
      obj.hs( obj.x == obj.bl ) = 0;

      % make variables at upper bound non basic
      obj.hs( obj.x == obj.bu ) = 1;

      % set fixed variables
      obj.hs( obj.bl == obj.bu ) = 4;

      % first m variables are basic
      obj.hs(1:obj.m) = 3;

    end

    function basis_main(obj)
      %basis_main  entry point to basis factorization
      %
      % Method calls:
      %   basis_activate
      %   basis_update
      %
      % On entry:
      %   stpbnd = {-1,1} if bound is hit, 0 if no bound is hit
      %   r_ix = index of limiting bound
      %   [other properties required for basis_update and basis_activate]
      %
      % On exit:
      %   upflg = true if factorization update is performed
      %   reflg = true if refactoriation is performed
      %   bsflg = true if bs factorization is performed
      %   brflg = true if br factorization is performed
      %   [other properties are modified by methods called]
      %

      % set all factorization flags to false
      obj.upflg = false;
      obj.reflg = false;
      obj.bsflg = false;
      obj.brflg = false;

      % reset b_ix
      % we don't yet know if a variable needs to move to the basic set
      obj.b_ix = 0;

      if obj.stpbnd ~= 0 && obj.hs(obj.r_ix) == 2
        % variable r_ix is superbasic and reaches bound
        % no change of basis factorization needed
        % just need to set r_ix as nonbasic in bs
        obj.basis_activate();
      elseif obj.stpbnd ~=0 && obj.hs(obj.r_ix) == 3
        % variable r_ix is basic and reaches bound
        % update to basis factorization needed
        obj.basis_update();
      end

      % note if obj.stpbnd == 0, there is no limiting variable.  Nothing
      % needs to be done.

    end

    function basis_update(obj)
      %update_basis  Update the basis factorization after a bound is hit
      %
      % This method is called when a basic variable reaches a bound.  The
      % basis factorization must be updated so that:
      %   - the column corresponding to the limiting variable is removed
      %   - a new column corresponding to a superbasic variable is added
      % If the update leads to a singular basis, a basis repair procedure
      % is envoked.
      %
      % (TODO) replace direct calls to bfac
      % (TODO) document method of selecting column from S
      %
      % Method calls:
      %   bfac.solveAt (no side effect)
      %   basis_activate
      %   fac_main
      %   bfac.repcol (changes basis factorization)
      %
      % On entry:
      %   !! many methods called
      %
      % On exit:
      %   !! many methods called
      %

      if sum(obj.hs == 2) == 1
        % there is only 1 superbasic variable, it must be used
        obj.b_ix = find(obj.hs == 2);
      elseif sum(obj.hs == 2) > 1
        % there is more than 1 superbasic variable, must choose one
        e = zeros(obj.m,1);
        e(obj.bcol == obj.r_ix) = 1;
        pi = obj.bfac.solveAt(e);
        w = zeros(obj.m+obj.n,1);
        w(obj.hs == 2) = abs(obj.A(:,obj.hs == 2)'*pi);
        [junk obj.b_ix] = max(w);
      end

      % s_ix is now the index of the column we wish to add to the basis
      obj.hs(obj.b_ix) = 3;

      % update the hs vector to make r_ix nonbasic
      obj.basis_activate();

      % perform operations on the basis factorization
      if obj.facucnt >= obj.opt.facfrq
        % refactorization frequency has been met
        % perform full refactorization
        obj.fac_main();
      elseif obj.b_ix ~= obj.r_ix
        % x(s_ix) becomes basic and x(r_ix) becomes nonbasic
        obj.upflg = true;

        % get index of column to replace in B
        jrep = find(obj.bcol == obj.r_ix);
        obj.bcol(jrep) = obj.b_ix;
        inform = obj.bfac.repcol(obj.A(:,obj.b_ix),jrep);
        obj.facucnt = obj.facucnt + 1;

        if inform == 8
          error('arcopt_nm_lc:basis_update','jrep is not between 1 and m.')
        end

        if inform ~= 0
          % basis update seemed to have problem
          % perform refactorization to fix the problem
          obj.fac_main()
        end
      end

    end

    function basis_activate(obj)
      %basis_activate  set basis information for limiting bound
      %
      % Sets basis information vector hs when variable r_ix hits a bound
      % and becomes nonbasic.
      %
      % On entry:
      %   r_ix = index of limiting bound, to become nonbasic
      %   stpbnd = -1 if bl is hit, 1 if bu is hit
      %   hs = basis information vector
      %
      % On exit:
      %   hs(r_ix) = set to 0 if r_ix nonbasic at bl
      %              set to 1 if r_ix nonbasic at bu
      %

      % variable obj.r_ix becomes nonbasic
      if obj.stpbnd == 1
        obj.hs(obj.r_ix) = 1;
      elseif obj.stpbnd == -1
        obj.hs(obj.r_ix) = 0;
      else
        error('arcopt_nm_lc:basis_activate','incorrect stpbnd.');
      end

    end

    %% generic methods
    function restore(obj,tol)
      %restore  move infeasible superbasics and nonbasics to their bounds
      %
      % This method will move all infeasible superbasic and nonbasic
      % variables to their bounds.  Infeasible superbasic variables will
      % become nonbasic at their bounds.  The basis information vector will
      % be set accordingly.
      %
      % On entry:
      %   m, n = problem size
      %   x = variable vector
      %   bl, bu = working bounds on the variables
      %   hs = basis information vector
      %
      % On exit:
      %   x = superbasics and nonbasics are now feasible wrt bl and bu
      %   hs = infeasible superbasics have become nonbasic
      %

      if nargin < 2
        tol = 0;
      end

      for i = 1:obj.m+obj.n
        if obj.hs(i) ~= 3

          % if variable is out of bounds, make feasible and nonbasic
          if obj.x(i) < obj.bl(i)-tol
            obj.x(i) = obj.bl(i)-tol;
            obj.hs(i) = 0;
          elseif obj.x(i) > obj.bu(i)+tol
            obj.x(i) = obj.bu(i)+tol;
            obj.hs(i) = 1;
          end

          % if variable is fixed, specify in hs
          if obj.bl(i) == obj.bu(i)
            obj.hs(i) = 4;
          end

        end
      end

    end

    function comp_basics(obj)
      %comp_basics  compute value of basic variables
      %
      % Give set values for all superbasics and nonbasics, compute the
      % value of all basic variables.  This requires a full rank,
      % factorized basis.
      %
      % On entry:
      %   x = variable vector with appropriate superbasics and nonbasics
      %   hs = basis information vector
      %   bfac = full rank, factorized basis object
      %   bcol = basis column map
      %   wm = workspace of size m
      %
      % On exit:
      %   x = basic variables now set
      %

      % (TODO) replace call to bfac.solveA
      obj.wm = -obj.bfac.solveA(obj.A(:,obj.hs~=3)*obj.x(obj.hs~=3));
      obj.x(obj.bcol) = obj.wm;
    end

    function basic_feas_chk(obj,tol)

      if nargin < 2
        tol = obj.exptol;
      end

      obj.bl_lx = ( obj.x < obj.bl_in-tol );
      obj.bu_lx = ( obj.x > obj.bu_in+tol );
      obj.infcnt = sum(obj.bl_lx) + sum(obj.bu_lx);
      obj.infsum = sum(obj.bl_in(obj.bl_lx)-obj.x(obj.bl_lx)) + ...
        sum(obj.x(obj.bu_lx)-obj.bu_in(obj.bu_lx));

    end

    function comp_y(obj)
      %comp_y  compute the vector of multipliers, y
      obj.y = obj.bfac.solveAt(obj.g(obj.bcol));
    end

    function comp_z(obj)
      %comp_z  compute the residual gradient vector, z
      obj.z(obj.hs == 3) = 0;
      obj.z(obj.hs ~= 3) = obj.g(obj.hs ~= 3) - obj.A(:,obj.hs ~= 3)'*obj.y;
    end

    function comp_dx_stpmax(obj)
      %comp_dx_stpmax compute max step size along dx vector
      %
      % On entry:
      %   dx = search direction
      %   x = current point
      %   bl = lower bound
      %   bu = upper bound
      %   opt.expsml = zero tolerance for expand
      %   exptol = working feasibility tolerance
      %   expinc = expand increase value
      %   opt.stpmax = max step size possible
      %   opt.infval = infinity value
      %
      % On exit:
      %   stpmax = max step size
      %   r_ix = index of limiting bound, 0 if no bound is hit
      %   stpbnd = -1 if bl hit, 1 if bu hit, 0 otherwise
      %   dgnflg = degenerate step if dgnflg == 1
      %

      idx = find(obj.dx ~= 0);
      [obj.stpmax obj.r_ix obj.stpbnd obj.dgnflg] = expand(idx,...
        obj.x, ...
        obj.dx, ...
        obj.bl, ...
        obj.bu, ...
        obj.opt.expsml, ...
        obj.exptol, ...
        obj.expinc, ...
        obj.opt.infval);

      % must check against stpmax parameter
      if obj.stpmax > obj.opt.stpmax
        % stpmax is too large, a bound will not be hit
        obj.stpmax = obj.opt.stpmax;
        obj.stpbnd = 0;
        obj.r_ix = 0;
      end

    end

    %% phase 1 methods
    function phase1_price(obj)
      %phase1_price  price nonbasic variables for phase1
      %
      % This method chooses the variable to leave the nonbasic set.
      %
      % On entry assumptions
      %   obj.z is the reduced gradient vector
      %   obj.hs is the basis information vector
      %
      % On finish
      %   obj.s_ix = index of variable to leave nonbasic set
      %   obj.sigma = direction to move variable obj.s_ix
      %   obj.hs(obj.s_ix) = 2, variable obj.s_ix is made superbasic
      %

      z_check = obj.z;
      z_check(obj.hs == 3) = 0;
      z_check(obj.hs == 0) = min(0,z_check(obj.hs == 0));
      z_check(obj.hs == 1) = min(0,-z_check(obj.hs == 1));

      % ignore fixed variables
      z_check(obj.hs == 4) = 0;

      [max_z obj.s_ix] = max(abs(z_check));
      obj.hs(obj.s_ix) = 2;
      obj.sigma = -sign(obj.z(obj.s_ix));
    end

    function phase1_dx(obj)
      %phase1_dx  compute search direction dx for phase1
      %
      % On entry assumptions
      %   obj.dx = vector with length m+n
      %   obj.hs = basis information vector
      %   obj.wm = workspace vector of size m
      %   obj.sigma = direction to move variable obj.s_ix
      %   obj.A = constraint matrix
      %   obj.bfac = lusol object with factorized basis matrix
      %   obj.bcol = map between basis matrix and columns of A
      %
      % On finish
      %   obj.dx = the search direction
      %

      % zero out nonbasic variables
      obj.dx(obj.hs ~= 3) = 0;

      % compute change for basic variables
      obj.wm = -obj.sigma*(obj.bfac.solveA(obj.A(:,obj.s_ix)));
      obj.dx(obj.bcol) = obj.wm;

      % set the direction for the free variable
      obj.dx(obj.s_ix) = obj.sigma;
    end

    function phase1_stp(obj)
      %phase1_stp  compute step size for phase 1

      % first do standard step size calculation
      obj.comp_dx_stpmax();
      obj.stp = obj.stpmax;

      if obj.dgnflg
        % step is degenerate, update counter and set stp
        obj.dgn1cnt = obj.dgn1cnt + 1;
        obj.stp = obj.stpmax;
      else
        % step is not degenerate,

        % find step size that removes as many infeasibilities as possible, but
        % does not go any further
        stptmp = 0;
        stpfeas = 0;
        for i = 1:obj.n+obj.m
          stptmp = 0;
          if obj.bl_lx(i) && obj.dx(i) ~= 0.0
            stptmp = (obj.bl_in(i)-obj.x(i))/obj.dx(i);
          elseif obj.bu_lx(i) && obj.dx(i) ~= 0.0
            stptmp = (obj.bu_in(i)-obj.x(i))/obj.dx(i);
          else
            % do nothing
          end
          if stptmp > stpfeas
            stpfeas = stptmp;
          end
        end

        if stpfeas < obj.stpmax
          % stpfeas is less than limiting step

          % no constraint becomes active
          obj.r_ix = 0;
          obj.stpbnd = 0;

          % set the step size
          obj.stp = stpfeas;

        end

      end

    end

    function phase1flg = phase1(obj)
      %phase1  carry out phase 1 procedure
      %
      % This method attempts to find a feasible point to the system of
      % inequalities:
      %
      %   bl <= x <= bu
      %   A(:,1:n)*x(1:n) == x(n+1:n+m)
      %
      %   (the last m elements of x are slack variables)
      %
      % It uses the simplex method to minimize the sum of infeasibilities.
      % The linear objective will change each iteration.
      %
      % On entry:
      %   This method uses most of the object properties.
      %   Most importantly this method expects consistent:
      %   x = current iterate
      %   hs = basis information vector
      %   bfac = basis factorization object
      %   bcol = basis column map
      %
      % On exit:
      %   (x,hs) = feasible point or optimal infeasible point
      %
      % Output:
      %   phase1flg = phase 1 status flag
      %
      %
      % Phase 1 status flags:
      %   phase1flg = 101 if feasible point found
      %   phase1flg = 111 if optimal point found, problem infeasible
      %   phase1flg = 131 if itermax is reached
      %

      obj.mprint('iter','\n');
      obj.mprint('iter','Phase 1 initiated.\n');

      obj.bl = obj.bl_in;
      obj.bu = obj.bu_in;

      % compute value of basic variables
      obj.comp_basics();

      % initialize the print counter
      obj.prntcnt = 1;
      phase1flg = 100;
      while(1)

        % debugging code
        %if obj.itercnt == 156
        %  keyboard
        %end

        % check iterations
        if obj.itercnt > obj.opt.itermax
          phase1flg = 131;
          break;
        end

        % expand
        obj.exp_main();

        % check feasibility of x
        obj.basic_feas_chk();

        if obj.infcnt == 0;
          % x is feasible, phase 1 is complete
          phase1flg = 101;
          break;
        end

        % set objective and bounds for phase 1
        obj.g = obj.bu_lx - obj.bl_lx;
        obj.bl = obj.bl_in;
        obj.bu = obj.bu_in;
        obj.bl(obj.bl_lx) = -obj.opt.infval;
        obj.bu(obj.bu_lx) = obj.opt.infval;
        %obj.bl(obj.bl_lx) = obj.x(obj.bl_lx);
        %obj.bu(obj.bu_lx) = obj.x(obj.bu_lx);

        % set basis for fixed nonbasic variables
        obj.hs(obj.hs ~= 3 & (obj.bl == obj.bu)) = 4;

        % compute the multipliers, y
        obj.comp_y();

        % compute the reduced gradient, z
        obj.comp_z();

        % phase 1 optimality check
        tbl = max(min(obj.x-obj.bl,obj.z));
        tbu = max(min(obj.bu-obj.x,-obj.z));
        obj.termchk = max(tbl,tbu);
        if obj.termchk <= obj.opt.dtol;
          phase1flg = 111;
          break;
        end

        obj.phase1_price();
        obj.phase1_dx();
        obj.phase1_stp();

        % take step
        obj.x = obj.x + obj.stp*obj.dx;

        % run basis routines
        obj.basis_main();

        % increment counters
        obj.phase1cnt = obj.phase1cnt + 1;
        obj.itercnt = obj.itercnt + 1;

        % print information for this iteration
        obj.phase1_print();

      end

      % print final phase 1 information
      obj.mprint('iter','\n');
      obj.mprint('iter','Phase 1 complete: ');
      switch phase1flg
        case 111
          obj.mprint('iter','problem infeasible.\n');
        case 101
          obj.mprint('iter','feasible point found.\n');
        case 131
          obj.mprint('iter','maximum number of iterations reached.\n');
        otherwise
          obj.mprint('iter','unknown.\n');
      end

      % clean up
      obj.prntcnt = 1;
      obj.bl = obj.bl_in;
      obj.bu = obj.bu_in;
      obj.g(:) = 0;
      obj.y(:) = 0;
      obj.z(:) = 0;

    end

    function phase1_print(obj)

      if ~strcmp(obj.opt.print_level,'iter')
        % do nothing
        return;
      end

      % currently printing
      % item: iter infcnt infsum D   stp   bfac
      % head: iter infcnt infsum D   stp   bfac
      % hspc: %5s  %6s    %8s    %1s %7s   %4s
      % dspc: %5d  %6d    %8.1e  %1s %7.1e %4s

      if obj.prntcnt == 1
        % print header
        obj.mprint('iter','\n');
        obj.mprint('iter','%5s %6s %8s %1s %7s %4s\n',...
          'iter','infcnt','infsum','D','stp','bfac');
      end

      if obj.prntcnt == obj.opt.print_freq
        obj.prntcnt = 1;
      else
        obj.prntcnt = obj.prntcnt + 1;
      end

      % get d string
      ds = ' ';
      if obj.dgnflg, ds = 'd'; end

      % get bfac string
      fu = ' '; ff = ' '; fs = ' '; fr = ' ';
      if obj.upflg, fu = 'u'; end
      if obj.reflg, ff = 'f'; end
      if obj.bsflg, fs = 's'; end
      if obj.brflg, fr = 'r'; end

      obj.mprint('iter','%5d %6d %8.1e %1s %7.1e %4s\n',...
        obj.itercnt,obj.infcnt,obj.infsum,ds,obj.stp,[fu ff fs fr]);

    end

    %% multiply methods
    function Y = mulZt(obj,X)
      %mulZt  Y = Z'*X
      %
      % On entry:
      %   hs = basis information vector
      %   bcol = basis map
      %   bfac = basis factorization object
      %
      % Input:
      %   X = input matrix with m+n rows and arbitrary number of columns
      %
      % Output:
      %   Y = Z'*X, number of rows depends on size of nullspace
      %             will have same number of columns as X
      %

      Y = X(obj.hs==2,:) - obj.A(:,obj.hs==2)'*obj.bfac.solveAt(X(obj.bcol,:));

    end

    function Y = mulZ(obj,X)
      %mulZ  Y = Z*X
      %
      % On entry:
      %   n = number of variables
      %   m = number of constraints
      %   hs = basis information vector
      %   bcol = basis map
      %   bfac = basis factorization object
      %
      % Input:
      %   X = input matrix, number of rows depends on size of nullspace
      %                     arbitrary number of columns
      %
      % Output:
      %   Y = Z*X, has m+n rows and same number of columns as X
      %

      Xc = size(X,2);
      Y = zeros(obj.n+obj.m,Xc);
      Y(obj.hs == 2,:) = X;
      Y(obj.bcol,:) = -obj.bfac.solveA(obj.A(:,obj.hs==2)*X);

    end

    function Y = mulZtHZ(obj,X)
      %mulZtHZ  compute product with reduced Hessian, Y = Z'*H*Z*X
      %
      % On entry:
      %   n = number of variables
      %   m = number of constraints
      %   H = Hessian matrix
      %
      % Input:
      %   X = input matrix, number of rows depends on size of nullspace
      %                     can have arbitrary number of columns
      %
      % Output:
      %   Y = Z'*H*Z*X, number of rows depends on size of nullspace
      %                 same number of columns as X
      %

      Y1 = obj.mulZ(X);
      Y1(1:obj.n,:) = obj.H*Y1(1:obj.n,:);
      Y1(obj.n+1:end,:) = 0;
      Y = obj.mulZt(Y1);

    end

    function Y = mulZtHZ_shift(obj,X,d)
      %mulZtHZ_shift  compute product with reduced Hessian shifted by d*I
      %
      % This method computes: Y = (Z'*H*Z + d*I)*X
      %
      % Input:
      %   X = input vector, number of rows depends on size of nullspace
      %                     arbitrary number of columns
      %   d = value of the shift
      %
      % Output:
      %   Y = (Z'*H*Z + d*I)*X
      %

      Y = obj.mulZtHZ(X) + d*X;

    end

    %% phase 2 methods
    function [s ds] = arc_eval(obj,stp)
      %arc_eval parameterized arc evaluation for use with arc search
      %
      % Input:
      %   stp = trial step size
      %
      % On entry:
      %   arc = trust region arc object
      %   arc_ZQ = arc nullspace matrix
      %
      % Output:
      %   s = arc position at stp
      %   d = derivative of arc at stp
      %

      s = obj.arc_ZQ*obj.arc.sol_s(stp);
      ds = obj.arc_ZQ*obj.arc.sol_s(stp,1);

    end

    function comp_hess(obj)
      %comp_hess  evaluate and store Hessian at current iterate
      %
      % Evaluates the Hessian of the objective function.  Note that it stores H
      % as a n by n sparse matrix.  Thus it is not possible to take products
      % with a n+m vector (containing slacks).
      %
      % On entry:
      %   n = number of true variables
      %   x = current iterate
      %
      % On exit:
      %   H = sparse Hessian at x with size n by n
      %

      obj.H = obj.usrhess(obj.x(1:obj.n));
    end

    function comp_eigs(obj)
      %comp_eigs  compute direction of negative curvature and properties
      %
      % This method uses Matlab's eigs function to compute the eigenvector
      % corresponding to the smallest eigenvalue.  This is used as the
      % direction of negative curvature for other parts of the algorithm.
      %
      % If there is a negative eigenvalue, then dnc is set to the the
      % eigenvector corresponding to the smallest eigenvalue.  If all
      % eigenvalues (according to eigs) are non-negative, then dnc = 0.
      %
      % Method calls:
      %   mulZtHZ
      %
      % On entry:
      %   hs = basis information vector
      %   opt.eigs_p = default number of storage vectors for eigs
      %   - Hessian and nullspace must be ready for computations
      %
      % On exit:
      %   evec = eigenvector
      %   eval = eigenvalue
      %   eigsflg = status flag from eigs, returns 0 if eigenvalue converged
      %   evalflg = true if negative eigenvalue is present
      %             false otherwise
      %

      % compute number of slack variables
      ns = sum(obj.hs == 2);

      % compute eigenvector corresponding to smallest eigenvalue
      if ns == 0
        % there are no superbasic variables at this point
        obj.evec = [];
        obj.eval = 0;
        obj.eigsflg = 0;
      elseif ns == 1
        % there is only one super basic variable
        % nothing to compute, just get entry of matrix
        obj.evec = 1;
        obj.eval = obj.mulZtHZ(1);
        obj.eigsflg = 0;
      else
        % there is more than one super basic variable

        % number of storage vectors for eigs
        obj.opt_eigs.p = min(ns,obj.opt.eigs_p);

        % starting vector for eigs
        obj.opt_eigs.v0 = ones(ns,1);

        % run eigs
        [obj.evec obj.eval obj.eigsflg] = eigs(@obj.mulZtHZ, ...
          ns, ...
          1, ...
          'sa', ...
          obj.opt_eigs);

        % if eigsflg = 0, then the eigenvalue converged
      end

      % set evalflg
      if ns > 0 && obj.eval < -obj.opt.ctol
        obj.evalflg = true;
      else
        obj.evalflg = false;
      end
      
    end

    function phase2_price(obj)
      %phase2_price determine if a nonbasic variable can be made superbasic
      %
      % Nothing is done if a negative eigenvalue is present.
      %
      % On entry:
      %   evalflg = true if there is a negative eigenvalue
      %   hs = basis information vector
      %   z = residual gradient vector
      %
      % On exit:
      %   s_ix = index of new superbasic variable, 0 if non selected
      %   hs = updated basis information vector
      %

      obj.s_ix = 0;

      % do nothing if there is a negative eigenvalue
      if obj.evalflg
        return;
      end
      
      % cycle_gaurd check
      if sum(obj.hs == 2) > 0 && ...
          max(abs(obj.z(obj.hs == 2))) >= obj.opt.dtol && ...
          obj.opt.cycle_guard && ...
          obj.stpbnd
        return;
      end
      
      % standard price
      if sum(obj.hs == 2) == 0 || ...
          max(abs(obj.z(obj.hs == 2))) < obj.opt.dtol

        z_check = obj.z;

        % ignore basic, superbasic, and fixed variables
        z_check(obj.hs == 2) = 0;
        z_check(obj.hs == 3) = 0;
        z_check(obj.hs == 4) = 0;

        z_check(obj.hs == 0) = min(0,z_check(obj.hs == 0));
        z_check(obj.hs == 1) = min(0,-z_check(obj.hs == 1));

        [max_z obj.s_ix] = max(abs(z_check));
        obj.hs(obj.s_ix) = 2;

      end

    end

    function comp_rg(obj)
      %comp_rg  compute reduced gradient at current iterate
      %
      % On entry:
      %   g = gradient
      %
      % On exit:
      %   rg = reduced gradient
      %

      obj.rg = obj.mulZt(obj.g);

    end

    function comp_dnc(obj)
      %comp_dnc  compute direction of negative curvature and properties
      %
      % The sign of the dnc is chosen such that rg'*dnc >= 0.
      %
      % On entry:
      %   hs = basis information vector
      %   rg = reduced gradient
      %   evalflg = negative eigenvalue flag, true if negative eigenvalue
      %             is present
      %   eval = eigenvalue
      %   evec = eigenvector
      %   dnccnt = dnc counter
      %
      % On exit:
      %   dnc = direction of negative curvature, zero vector if all eigenvalues
      %         are non-negative
      %   dncflg = true if dnc exists, false if not
      %   dnccnt = dnc counter, possibly incremented
      %

      % compute number of slack variables
      ns = sum(obj.hs == 2);

      % set dnc properties
      if obj.evalflg
        % dnc exists
        obj.dnc = obj.evec;
        obj.dncflg = true;
        obj.dnccnt = obj.dnccnt + 1;

        % set appropriate direction for dnc
        if obj.rg'*obj.dnc < 0
          obj.dnc = -obj.dnc;
        end

      else
        % dnc does not exist
        obj.dnc = zeros(ns,1);
        obj.dncflg = false;
      end

    end

    function phase2_dx(obj)
      %phase2_dx  compute steepest descent direction for phase 2
      %
      % Computes dx = -Z*Z'*g, which is the steepest descent direction on
      % on the current nullspace.
      %
      % Method calls:
      %   mulZtHZ
      %
      % On entry:
      %   rg = reduced gradient vector
      %   dncflg = true if dnc exists
      %   dnc = direction of negative curvature
      %   opt.dnc_mtol = tolerance of dnc modification
      %   dncmcnt = dnc modification counter
      %
      % On exit:
      %   dx = steepest descent direction, possibly modified
      %   rdx = reduced steepest descent direction, possibly modified
      %   dncmflg = true if dnc modification needed
      %   dncmcnt = dnc modification counter, possibly incremented
      %

      % by default the direction is just steepest descent
      obj.rdx = -obj.rg;

      % set the dnc modification flag to false initially
      obj.dncmflg = false;

      if obj.dncflg
        % check condition for dnc modification
        v1 = obj.mulZtHZ(obj.dnc);
        rg_check = obj.rg'*v1;
        dnc_check = obj.dnc'*v1;
        if abs(rg_check) <= abs(dnc_check)*obj.opt.dnc_mtol
          % dnc modification required
          obj.rdx = -(obj.rg + abs(obj.eval)*obj.dnc);
          %obj.rdx = -(obj.rg + obj.dnc);
          obj.dncmflg = true;
          obj.dncmcnt = obj.dncmcnt + 1;
        end
      end

      obj.dx = obj.mulZ(obj.rdx);

    end

    function comp_lin(obj)
      %comp_lin  compute modified newton search direction
      %
      % Solves
      %
      %    (Z'HZ + delta*I)*p = -rg
      %
      % with
      %
      %    delta = max(-eval+lin_delta,lin_delta).
      %
      % lin_delta is a user specified regularization parameter.
      %
      % Method calls:
      %   mulZtHZ_shift
      %
      % On entry:
      %   obj.eval = minimum eigenvalue of reduced hessian
      %   obj.rg = reduced gradient vector
      %   obj.opt.lin_delta = lin shift parameter
      %   obj.opt.lin_tol = lin convergence tolerance
      %   obj.opt.lin_maxit = maximum number of lin iterations
      %
      % On exit:
      %   p = modified newton direction
      %   linflg = lin convergence flag
      %   linres = lin relative residual
      %   lincnt = lin iteration count
      %

      % compute delta perturbation
      delta = max(-obj.eval+obj.opt.lin_delta,obj.opt.lin_delta);

      % compute maxit
      if isempty(obj.opt.lin_maxit)
        maxit = 50*sum(obj.hs == 2);
      else
        maxit = obj.opt.lin_maxit;
      end

      % solve shifted linear system with pcg
      %[obj.p obj.linflg obj.linres obj.lincnt] = ...
      %  pcg(@(x) obj.mulZtHZ_shift(x,delta),-obj.rg, ...
      %  obj.opt.lin_tol,maxit);

      % solve shifted linear system with minres
      [obj.p obj.linflg obj.linres obj.lincnt] = ...
        minres(@(x) obj.mulZtHZ_shift(x,delta),-obj.rg, ...
        obj.opt.lin_tol,maxit);

      % solve linear system with minres
      %[obj.p obj.linflg obj.linres obj.lincnt] = ...
      %  minres(@(x) obj.mulZtHZ(x),-obj.rg, ...
      %  obj.opt.lin_tol,maxit);

      % pcg flags:
      % 0) pcg converged to the desired tolerance tol within maxit iterations.
      % 1) pcg iterated maxit times but did not converge.
      % 2) Preconditioner M was ill-conditioned.
      % 3) pcg stagnated. (Two consecutive iterates were the same.)
      % 4) One of the scalar quantities calculated during pcg became too small or too
      %    large to continue computing.

    end

    function comp_arc(obj)
      %comp_arc  construct the arc object
      %
      % This method constructs the arc object.  The first question is the size
      % of the subspace.  If a direction of negative curvature exists, then the
      % arc subspace is chosen to be [rg p dnc], where rg is the reduced
      % gradient, p is the modified newton direction, and dnc is direction of
      % negative curvature.  If no direction of negative curvature exists, the
      % subspace is chosen to be [rg p].
      %
      % The right hand side of the arc must be perturbed under a condition in
      % order to prove convergence to a point satisfying the second order
      % optimality conditions.  The condition is
      %
      %   |g'*H*dnc| <= arc_mtol*|dnc'*H*dnc|                         (1)
      %
      %     g = gradient
      %     H = Hessian
      %     arc_mtol = modification tolerance
      %     dnc = direction of negative curvature
      %     (note that g, H, and dnc may be in reduced space)
      %
      % If (1) holds, then the right hand side of the arc is set to g+d.
      %
      % External classes used:
      %   Arc
      %
      % Method calls:
      %   mulZ
      %   mulZtHZ
      %
      % On entry:
      %   rg = reduced gradient vector, Z'*g
      %   p = modified newton direction
      %   dnc = direction of negative curvature
      %   dncflg = true indicates existence of dnc, false indicates otherwise
      %   hs = basis information vector, used to determine |S|
      %   rdx = reduced steepest descent direction, possibly modified
      %
      % On exit:
      %   arc_dim = dimension of the arc subspace
      %   arc_Q = orthonormal arc subspace matrix
      %   arc_ZQ = arc nullspace matrix
      %   arc_H = arc matrix
      %   arc_g = arc right hand side
      %   arc = arc object
      %

      % compute number of slack variables
      ns = sum(obj.hs == 2);

      % construct the subspace matrix
      if ns == 1
        % only one superbasic variable
        % use reduced gradient as subspace
        obj.arc_dim = 1;
        D = obj.rg;
      elseif ns == 2 && obj.dncflg
        % two superbasic variables
        %use reduced gradient and dnc as subspace
        obj.arc_dim = 2;
        D = [obj.rg obj.dnc];
      elseif ns > 2 && obj.dncflg
        % more than two superbasic variables
        % direction of negative curvature exists
        % subspace is spanned by [rg p dnc]
        obj.arc_dim = 3;
        D = [obj.rg obj.p obj.dnc];
      else
        % two or more superbasic variables
        % direction of negative curvature does not exist
        % subspace is spanned by [rg p]
        obj.arc_dim = 2;
        D = [obj.rg obj.p];
      end

      % compute orthogonal subspace matrix
      [obj.arc_Q R] = qr(D,0);

      % compute arc nullspace matrix
      obj.arc_ZQ = obj.mulZ(obj.arc_Q);

      % construct arc matrix
      obj.arc_H = obj.arc_Q'*obj.mulZtHZ(obj.arc_Q);

      % construct arc rhs
      obj.arc_g = -obj.arc_Q'*obj.rdx;

      % create arc
      % (TODO) reuse same arc object
      obj.arc = trarc(obj.arc_H,obj.arc_g);

    end

    function comp_stpinit(obj)
      %comp_stpinit  compute initial step size
      %
      % Constructs the initial step size for the s-parameterization of the
      % arc.
      %
      % If the Hessian is positive definite, the initial step size is:
      %   s0 = 1/min(eig(H))
      %
      % If the Hessian is positive semi-definite, the initial step size is:
      %   s0 = 1
      % (TODO) this is probably not a good choice.
      %
      % If the Hessian is indefinite, the initial step size is:
      %   s0 = -1/min(eig(H))
      % This choice was used in both [Behrman1998] and [Gatto2000].
      % (TODO) a better selection can be made here as well.
      %
      % On entry:
      %   arc = arc object
      %   arc.v_min = min eigenvalue for arc matrix, or min eigenvalue of Hessian
      %   opt.arc_vtol = tolerance for selecting trial step size
      %
      % On exit:
      %   obj.stpinit = initial step size
      %

      v_min = obj.arc.v_min;

      if v_min >= obj.opt.arc_vtol
        % arc is in positive definite subspace
        obj.stpinit = 1/v_min;
      elseif v_min <= -obj.opt.arc_vtol
        % arc is in an indefinite subspace
        obj.stpinit = -1/v_min;
      else
        % arc is in positive semi-definite subspace
        obj.stpinit = 1;
      end

    end

    function comp_arc_stpmax(obj)
      %comp_arc_stpmax  compute the maximum step size along arc
      %
      % On entry:
      %   arc = trust region arc
      %   x = current iterate
      %   arc_ZQ = arc nullspace matrix
      %   bl = lower bound vector
      %   bu = upper bound vector
      %   opt.stpmax = max value of stp
      %   exptol = feasibility tolerance from expand
      %   opt.infval = infinity value
      %
      % On exit:
      %   stpmax = max step size
      %   r_ix = index of limiting bound, 0 if step is unrestricted
      %   stpbnd = -1 if lower bound hit, 1 if upper bound hit, 0 if no bound hit
      %

      [obj.stpmax obj.r_ix obj.stpbnd] = trarc_exptst(obj.arc, ...
        obj.x, ...
        obj.arc_ZQ, ...
        obj.bl, ...
        obj.bu, ...
        obj.opt.stpmax, ...
        obj.exptol, ...
        obj.opt.infval);

    end

    function comp_srch_d2val(obj)
      %comp_srch_d2val  compute inital value for second derivative of search function
      %
      % Compute the srch_d2val parameter, which is need for the arc search.
      % Let's say the objective function is f(x) and the search arc is w(s).
      % The search function is:
      %
      %   phi(s) = f(x + w(s))
      %
      % The initial values of the first and second derivatives are:
      %
      %   phi'(0) = g'*w'(0)
      %   phi''(0) = w'(0)'*H*w'(0) + g'*w''(0)
      %
      % This method computes phi''(0) and stores it in srch_d2val.
      %
      % On entry:
      %   - the arc is fully initialized
      %   arc = trust region arc object
      %   arc_ZQ = arc nullspace matrix
      %   g = gradient of objective function
      %   H = Hessian of objective function
      %   n = number of true variables
      %
      % On exit:
      %   srch_d2val = initial value for second derivative of search function
      %

      % compute initial first derivative of arc
      [s0 ds0] = obj.arc_eval(0);

      % compute initial second derivative of arc
      d2s0 = obj.arc_ZQ*obj.arc.sol_s(0,2);

      % compute initial value for second derivative of search function
      obj.srch_d2val = obj.g'*d2s0 + ds0(1:obj.n)'*(obj.H*ds0(1:obj.n));

    end

    function phase2_srch(obj)
      %phase2_srch  arc search for phase 2

      % compute maximum function evaluations for this call to arc search
      fevmax = min(obj.opt.fevas,obj.opt.fevmax - obj.fevcnt);

      % perform arc search
      [obj.x obj.f obj.g obj.stp obj.srch_flg obj.srch_fev] = mtasrch(...
        @obj.func_wrap, ...
        obj.x, ...
        obj.f, ...
        obj.g, ...
        @obj.arc_eval, ...
        obj.stpinit, ...
        obj.srch_d2val, ...
        obj.opt.ftol, ...
        obj.opt.gtol, ...
        obj.opt.xtol, ...
        obj.opt.stpmin, ...
        obj.stpmax, ...
        fevmax);

      % increment function evaluation counter
      obj.fevcnt = obj.fevcnt + obj.srch_fev;

      if obj.stpbnd == 0 || ...
          ( obj.srch_flg ~= 2 && obj.srch_flg ~= 5 && obj.srch_flg ~= 6 )
        % a bound was not hit
        obj.stpbnd = 0;
        obj.r_ix = 0;
      end

    end

    function termflg = phase2_term(obj)
      %phase2_term termination checks for phase 2
      %
      % Check termination conditions for phase 2.
      %
      % Termination flags:
      %  0 = no termination, algorithm should continue
      %  1 = optimality conditions satisfied
      %  2 = reached maximum number of function evaluations
      %  3 = reached maximum number of iterations
      %  4 = arc search reported a failure
      %

      % no termination
      termflg = 0;

      % check the number of function evaluations
      if obj.fevcnt >= obj.opt.fevmax
        termflg = 2;
      end

      % check iterations
      if obj.itercnt > obj.opt.itermax
        termflg = 3;
      end

      % check arc search termination flag
      %   srch_flg = 1 if stp satisfies the descent and curvature condition
      %   srch_flg = 2 if interval size is less than xtol
      %   srch_flg = 3 if algorithm has exceeded maxfev
      %   srch_flg = 4 if stpmin > 0 and stp == stpmin
      %   srch_flg = 5 if stp == stpmax
      %   srch_flg = 6 if stp == stpmax & strong wolfe conditions met
      %   srch_flg = 7 if rounding errors prevent progress
      if ismember(obj.srch_flg,[2 3 4 7])
        termflg = 4;
      end

      % phase 2 optimality check
      tbl = max(min(obj.x-obj.bl,obj.z));
      tbu = max(min(obj.bu-obj.x,-obj.z));
      obj.termchk = max(tbl,tbu);
      if obj.termchk <= obj.opt.dtol && obj.eval >= -obj.opt.ctol;
        termflg = 1;
      end

    end

    function phase2flg = phase2(obj)
      %phase2 carry out phase 2 of algorithm
      %
      % phase 2 exit flags
      %   phase2flg = 201 if optimal
      %               231 if itermax is hit
      %               232 if fevmax is hit

      obj.mprint('iter','\n');
      obj.mprint('iter','Phase 2 initiated.\n');

      % reset basis change indices
      obj.s_ix = 0;
      obj.r_ix = 0;
      obj.b_ix = 0;
      obj.stpbnd = 0;

      % reset srch_flg
      obj.srch_flg = 0;

      % evaluate function
      obj.g = zeros(obj.n+obj.m,1);
      [obj.f obj.g] = obj.func_wrap(obj.x);
      obj.fevcnt = 1;

      obj.prntcnt = 1;
      phase2flg = 0;
      while(1)

        %if obj.itercnt == 100
        %  keyboard
        %end

        %fprintf('iter=%d hs(8)=%d\n',obj.itercnt,obj.hs(8))

        % save old values
        obj.x_old = obj.x;
        obj.f_old = obj.f;
        
        % expand
        obj.exp_main();

        % compute multiplier vector, y
        obj.comp_y();

        % compute residual gradient vector, z
        obj.comp_z();

        % evaluate hessian
        obj.comp_hess();

        % compute direction of negative curvature
        obj.comp_eigs();

        % check termination conditions
        phase2flg = phase2_term(obj);

        % terminate if phase2flg ~= 0
        if phase2flg
          break;
        end

        % decide whether to delete constraint
        obj.phase2_price();

        % compute reduced gradient vector, rg
        obj.comp_rg();

        % construct dnc
        obj.comp_dnc();

        % compute the steepest descent direction in nullspace
        obj.phase2_dx();

        % compute the maximum step size along dx
        obj.comp_dx_stpmax();

        % increment phase 2 degeneracy counter if needed
        if obj.dgnflg, obj.dgn2cnt = obj.dgn2cnt + 1; end

        if obj.dgnflg
          % at degenerate point, take step and evaluate
          obj.stp = obj.stpmax;
          obj.x = obj.x + obj.stp*obj.dx;

          obj.g(1:obj.n+obj.m) = 0;
          [obj.f obj.g(1:obj.n)] = obj.usrfun(obj.x(1:obj.n));
          obj.fevcnt = obj.fevcnt + 1;

          % set some flags for the methods that won't be evaluated
          obj.linflg = 0;
          obj.lincnt = 0;
          obj.linres = 0;

        else
          % non-degenerate point, proceed with algorithm

          % compute the modified newton direction
          obj.comp_lin();

          % construct the arc object
          obj.comp_arc();

          % compute the initial step length
          obj.comp_stpinit();

          % compute the max step size along arc with respect to constraints
          obj.comp_arc_stpmax();

          % compute initial value of second derivative of search function
          obj.comp_srch_d2val();

          % carry out arc search
          obj.phase2_srch();

        end

        % update basis and factorization if needed
        obj.basis_main();

        % increment counters
        obj.itercnt = obj.itercnt + 1;
        obj.phase2cnt = obj.phase2cnt + 1;

        % print information for this iteration
        obj.phase2_print();

      end

      % print final phase 2 information
      obj.mprint('iter','\n');
      obj.mprint('iter','Phase 2 complete: ');
      switch phase2flg
        case 1
          obj.mprint('iter','optimality conditions satisfied.\n');
        case 2
          obj.mprint('iter','maximum number of function evaluations reached.\n');
        case 3
          obj.mprint('iter','maximum number of iterations reached.\n');
        otherwise
          obj.mprint('iter','unknown.')
      end

    end

    function phase2_print(obj)

      if ~strcmp(obj.opt.print_level,'iter')
        % do nothing
        return;
      end

      % currently printing
      % item: iter fevcnt f     ||z|| D   eval  M   s   scnt sres  stp   bas bfac
      % head: iter fevcnt f     ||z|| D   eval  M   s   scnt sres  stp   bas bfac
      % hspc: %5s  %6s    %8s   %7s   %1s %8s   %1s %1s %4s  %7s   %7s   %3s %4s
      % dspc: %5d  %6d    %8.1e %7.1e %1s %8.1e %1s %1d %4d  %7.1e %7.1e %3s %4s

      % item: iter ||z|| D   eval  M   s   prd sres  stp   bas bfac fev df    dx    ns
      % head: iter ||z|| D   eval  M   s   prd sres  stp   bas bfac fev df    dx    ns
      % hspc: %4s  %5s   %1s %6s   %1s %1s %3s %5s   %5s   %3s %4s  %3s %6s   %5s   %s 
      % dspc: %4d  %5.0e %1s %6.0e %1s %1d %3d %5.0e %5.0e %3s %4s  %3d %6.0e %5.0e %d

      
      if obj.prntcnt == 1
        % print header
        obj.mprint('iter','\n');
        obj.mprint('iter','iter: %d, fev: %d, f: %6.0e, ||x||: %5.0e\n\n',obj.itercnt,obj.fevcnt,obj.f,norm(obj.x,'inf'));
        obj.mprint('iter','%4s %5s %1s %6s %1s %1s %3s %5s %5s %3s %4s %3s %6s %5s %s\n',...
          'iter','||z||','D','eval','M','s','prd','sres','stp','bas','bfac','fev','df','dx','ns');
      end

      if obj.prntcnt == obj.opt.print_freq
        obj.prntcnt = 1;
      else
        obj.prntcnt = obj.prntcnt + 1;
      end

      % handle eigenvalue
      if obj.dgnflg
        tmp_eval = 0;
      else
        tmp_eval = obj.eval;
      end

      % get m string
      ms = ' ';
      if obj.dncmflg, ms = 'm'; end

      % get d string
      ds = ' ';
      if obj.dgnflg, ds = 'd'; end

      % get bas string
      bs = ' '; bn = ' '; bb = ' ';
      if obj.s_ix, bs = 's'; end
      if obj.r_ix, bn = 'n'; end
      if obj.b_ix, bb = 'b'; end

      % get bfac string
      fu = ' '; ff = ' '; fs = ' '; fr = ' ';
      if obj.upflg, fu = 'u'; end
      if obj.reflg, ff = 'f'; end
      if obj.bsflg, fs = 's'; end
      if obj.brflg, fr = 'r'; end

      obj.mprint('iter','%4d %5.0e %1s %6.0e %1s %1d %3d %5.0e %5.0e %3s %4s %3d %6.0e %5.0e %d\n',...
        mod(obj.itercnt,9999),...
        obj.termchk,...
        ds, ...
        tmp_eval,...
        ms,...
        obj.linflg,...
        obj.lincnt,...
        obj.linres,...
        obj.stp,...
        [bs bn bb],...
        [fu ff fs fr],...
        obj.srch_fev,...
        obj.f-obj.f_old,...
        norm(obj.x-obj.x_old,'inf'),...
        sum(obj.hs == 2));

    end

    %% debugging methods
    function dbg_plot_farc(obj,stp)
      stp = stp(:)';
      n_stp = length(stp);
      [w dw] = obj.arc_eval(stp);
      f = arrayfun(@(i) obj.func_wrap(obj.x+w(:,i)),1:n_stp);
      plot(stp,f);
      keyboard
    end

    function dbg_plot_fgd(obj,stp)
      stp = stp(:)';
      [w dw] = obj.arc_eval(0);
      f = arrayfun(@(alpha) obj.func_wrap(obj.x-alpha*obj.g),stp);
      plot(stp,f);
    end

    function dbg_status(obj)
      fprintf('n = %d\n',obj.n);
      fprintf('m = %d\n',obj.m);
      fprintf('size(A) = [%d,%d]\n',size(obj.A));
      fprintf('sum(hs==-1) = %d\n',sum(obj.hs==-1));
      fprintf('sum(hs==0) = %d\n',sum(obj.hs==0));
      fprintf('sum(hs==1) = %d\n',sum(obj.hs==1));
      fprintf('sum(hs==2) = %d\n',sum(obj.hs==2));
      fprintf('sum(hs==3) = %d\n',sum(obj.hs==3));
      fprintf('sum(hs==4) = %d\n',sum(obj.hs==4));
    end

    function dbg_state(obj,hs_set)

      if nargin == 1
        hs_set = [-1 0 1 2 3 4];
      end

      fprintf('%4s %2s %8s %8s %8s %8s\n', ...
        'idx','hs','bl','x','bu','g');

      for i = 1:obj.n+obj.m
        if ismember(obj.hs(i),hs_set)
          fprintf('%4d %2d %8.1e %8.1e %8.1e %8.1e\n', ...
            i,obj.hs(i),obj.bl(i),obj.x(i),obj.bu(i),obj.g(i));
        end
      end

    end
    
    function dbg_basis_check(obj)
      %dbg_basic_check  basis consistency check
      
      ptol = obj.opt.ptol;
      
      fprintf('number of variables, n = %d\n',obj.n);
      fprintf('number of constraints, m = %d\n',obj.m);
      
      for i = [-1 0 1 2 3 4]
        fprintf('sum(hs == %d) = %d\n',i,sum(obj.hs==i));
      end
      
      check = true(obj.n+obj.m,1);
      
      for i = 1:(obj.n + obj.m)
        switch obj.hs(i)
          case -1 % nonbasic between bl and bu
            if obj.x(i) < obj.bl_in(i) - ptol || obj.x(i) > obj.bu_in(i) + ptol
              check(i) = false;
            end
          case 0 % nonbasic at bl
            if obj.x(i) > obj.bl_in(i) + ptol
              check(i) = false;
            end
          case 1 % nonbasic at bu
            if obj.x(i) < obj.bu_in(i) - ptol
              check(i) = false;
            end
          case 2 % superbasic
            if obj.x(i) < obj.bl_in(i) - ptol || obj.x(i) > obj.bu_in(i) + ptol
              check(i) = false;
            end
          case 3 % basic
            if obj.x(i) < obj.bl_in(i) - ptol || obj.x(i) > obj.bu_in(i) + ptol
              check(i) = false;
            end
          case 4 % nonbasic, fixed bl == bu
            if obj.bl_in(i) ~= obj.bu_in(i) || abs(obj.x(i)-obj.bl_in(i)) > ptol
              check(i) = false;
            end
          otherwise
            check(i) = false;
            fprintf('unknown basis indicator\n');
        end
      end
      
      incon = ~check;
      incon_idx = find(incon);
      
      fprintf('number of inconsistent variables = %d\n',sum(incon));
      fprintf('inconsisten variables: ');
      for i = 1:length(incon_idx)
        fprintf('%d ',incon_idx(i));
      end
      fprintf('\n');
      
    end

  end

end
