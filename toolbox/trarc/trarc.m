%trarc  trust region arc class
%
% This class defines trust region arc objects.
%
% The trust region subproblem is:
%  
%   min 1/2 * w*H*w + g'*w
%   s/t w'*w <= c
%
% The optimality conditions are:
%
%    (H+y*I)*w == -g
%   y*(w'*w-c) == 0
%         w'*w <= c
%            y >= 0
%      (H+y*I) >= 0 (matrix is positive semi-definite)
%
% Here, y is the lagrange multiplier corresponding to the constraint.
%
% As c is increased (y is decreased) the solution to the problem traces out an
% arc.  This class provides a way to evaluate points along the arc with an
% eigendecomposition of H and then reparameterizing with
%
%   s(y) = 1 / (v_min + y)
%
% Where v_min is the smallest eigenvalue of H.  Given an eigendecomposition of
% H = U*V*U' the arc, parameterized by y, is
%
%   w(y) = -U*qy(y,V)*U'*g.
%
% Here, 
%
%   qy(y,v_i) = 1 / (v_i + y).
%
% When we reparameterize with s we get
%
%   w(s) = -U*qs(s,V)*U'*g,
%   qs(s,v_i) = s / (s*(v_i-v_min) + 1).
%
% trarc objects are also able to compute intersections with linear constraints
% and find step lengths for solutions with a given norm.  See methods:
%
%   trarc.intersect - intersection with linear constraint
%   trarc.bound - get step length for solution with given norm
%
% Example:
%   n = 2;
%   H = randn(n,n);
%   H = 0.5*(H+H');
%   g = randn(n,1);
%   arc1 = trarc(H,g);
%   s = linspace(0,5,20);
%   w = arc1.sol_s(s);
%   plot(w(1,:),w(2,:))
%

classdef trarc < handle
  
  properties (SetAccess = private)
    
    % basic properties
    H = []; % Hessian
    g = []; % gradient
    U = []; % eigenvector matrix
    v = []; % eigenvalues
    n = 0; % problem size
    utg = []; % storage for U'*g
    v_min = []; % minimum eigenvalue
    v_min_ix = []; % index of minimum eigenvalue
    v_max = []; % maximum eigenvalue
    v_max_ix = []; % index of maximum eigenvalue
    v_hat = []; % v - v_min
    
    % intersection with linear constraint
    int_flg = false; % true if int_pset is initialized, false otherwise
    int_pset = []; % set of polynomials for linear constraint
    int_p = []; % final polynomial
    
    % intersection with the trust region bound
    bnd_flg = false; % true if bnd_pset is initialized, false otherwise
    bnd_pset = []; % set of polynomials for tr bound constraint
    bnd_p = []; % final polynomial
    
  end
  
  methods
    
    function obj = trarc(H_in,g_in)
      %trarc constructor for trarc class.
      %
      % Calls set method on input Hessian and gradient.
      %
      
      obj.set(H_in,g_in);
    end
    
    function set(obj,H_in,g_in)
      %set set trarc object data
      %
      % Sets the problem size (obj.n) and calls setter method for H and g.
      %
      
      obj.n = size(H_in,1);
      obj.set_H(H_in);
      obj.set_g(g_in);
      
    end
    
    function set_H(obj,H_in)
      %set_H set Hessian data
      %
      % Checks size and symmetry of H.  Computes eigendecomposition and clears
      % out polynomial data.
      %

      % check size of H
      [n1 n2] = size(H_in);
      if n1 ~= obj.n || n2 ~= obj.n
        error('trarc:set_H','H does not have size [n n].');
      end
      
      % check that H is symmetric
      sym_err = max(max(abs(H_in-H_in')));
      if sym_err > 1e-8*(norm(H_in,inf) + 1)
        error('trarc:set_H','H must be a symmetric matrix. Error in symmetry is %g.',sym_err);
      end
      
      % set problem data
      obj.H = 0.5*(H_in+H_in');
      
      % eigendecomposition
      [obj.U V] = eig(obj.H);
      obj.v = diag(V);
      [obj.v_min obj.v_min_ix] = min(obj.v);
      [obj.v_max obj.v_max_ix] = max(obj.v);
      obj.v_hat = obj.v - obj.v_min;
      
      % not prepared for computing intersections
      obj.int_flg = false;
      obj.int_pset = [];
      obj.int_p = [];
      
      % not prepared for trust region bound calculation
      obj.bnd_flg = false;
      obj.bnd_pset = [];
      obj.bnd_p = [];

    end
    
    function set_g(obj,g_in)
      %set_g set gradient data
      %
      % Note this method must be called after H is set, because it makes use of
      % the eigenvector matrix.
      %
      
      % orient g
      obj.g = g_in(:);

      % check length
      if length(obj.g) ~= obj.n
        error('trarc:set_g','g_in does not have length n.');
      end
      
      % utg
      obj.utg = obj.U'*obj.g;

    end
    
    function q = qs(obj,s,d)
      %qs  evaluate function that operates on eigenvalues
      %
      % Evaluates the parameterized function that operates on eigenvalues for
      % the trust region arc.  We have v_hat = v - v_min, where v are
      % eigenvales of H and v_min is min(v).  The expressions for the function
      % and derivates are:
      %
      %   q(s,i) = s / (s*v_hat_i + 1)
      %   q'(s,i) = 1 / (s*v_hat_i + 1)^2
      %   q''(s,i) = -2 * v_hat_i / (s*v_hat_i + 1)^3
      %
      % Here i indicates the index of the eigenvalue under consideration.
      %
      % These derivatives may be checked with the following matlab code:
      %
      %   syms s v v_min real
      %   v_s = solve('s = 1/(v_min+v)',v)
      %   f1 = 1/(v+v_s)
      %   f = s / (s*(v-v_min)+1)
      %   df = simplify(diff(f,s))
      %   d2f = simplify(diff(df,s))
      %
      % Input:
      %   s = vector of parameter values, has length m
      %   d = derivative level {0,2,1}
      %
      % Output:
      %   q = function or derivative value for all eigenvalues, has size n by m
      %       where n is size of the problem and m is the length of s
      %
      
      % optional input
      % return function value if d not specified
      if nargin < 3
        d = 0;
      end

      % orient and get size of input
      s = s(:)';
      m = length(s);
      
      % compute q
      switch d
        case 0
          % evaluate function:
          %   q(s,i) = s / (s*v_hat_i + 1)
          q = repmat(s,obj.n,1) ./ (obj.v_hat*s + 1);
        case 1
          % evaluate first derivative:
          %   q'(s,i) = 1 / (s*v_hat_i + 1)^2
          q = 1 ./ (obj.v_hat*s + 1).^2;
        case 2
          % evaluate second derivative:
          %   q''(s,i) = -2 * v_hat_i / (s*v_hat_i + 1)^3
          q = -2*repmat(obj.v_hat,1,m) ./ (obj.v_hat*s + 1).^3;
        otherwise
          error('trarc:qs','d not in {0,1,2}. Cannot evaluate higher derivatives.')
      end
      
    end
    
    function w = sol_s(obj,s,d)
      %sol_s  produce solution parameterized by s
      %
      % Computes the solution:
      %
      %   w = -U*qs(s,v)*U'*g
      %
      % Input:
      %   s = vector of parameter values, has length m
      %   d = derivative level {0,2,1}
      %
      % Output:
      %   w = solution or derivative values, has size n by m
      %       where n is size of the problem and m is the length of s
      %
      
      % optional input
      % return function value if d not specified
      if nargin < 3
        d = 0;
      end

      % orient input
      s = s(:)';
      m = length(s);
      
      % first get q values
      w = obj.qs(s,d);
      
      % compute solution
      w = -obj.U * (w .* repmat(obj.utg,1,m));

    end
    
    function y = mult(obj,s)
      %mult  return trust region multiplier for given values of s
      %
      % The trust region multiplier is:
      %
      %   y = 1/s - v_min
      %
      % where v_min is the minimum eigenvalue of H.
      %
      % Input:
      %   s = vector of parameter values, has length m
      %
      % Output:
      %   y = multiplier values for given s, has length m
      %
      
      % orient s
      s = s(:)';
      
      % compute multiplier values
      y = 1./s - obj.v_min;
      
    end
    
    function s_int = intersect(obj,a,b)
      %intersect  compute the intersection between tr arc and linear constraint
      %
      % Computes the intersection of the arc with a linear constraint:
      %
      %   a'*w(s) >= b.
      %
      % The method assumes that w(0) is feasible, or b < 0.  The method returns
      % the smallest s_int > 0 such that:
      %
      %   a'*w(s_int) == b.
      %
      % If there is no intersection, s_int is returned as empty.
      %
      % Input:
      %   a = constraint vector
      %   b = constratin rhs
      %
      % Output:
      %   s_int = step size of intersection, is empty if no intersection with
      %           s_int > 0.
      %
      
      % first check to see if polynomial data needs to be initialized
      % the initialization need only be done once for a given H
      if ~obj.int_flg
        % need to construct set of polynomials
        obj.int_pset = zeros(obj.n+1,obj.n+1);
        lx = true(1,obj.n);
        for i = 1:obj.n
          lx(i) = false;
          obj.int_pset(i,2:obj.n+1) = poly(-obj.v_hat(lx));
          lx(i) = true;
        end
        obj.int_pset(obj.n+1,:) = poly(-obj.v_hat);
        obj.int_p = zeros(1,obj.n+1);
        obj.int_flg = true;
      end
      
      % compute coefficients
      a = a(:);
      alpha = -(obj.U'*a) .* obj.utg;
      
      % reset the polynomial
      obj.int_p(:) = 0;

      % construct polynomial for this intersection
      for i = 1:obj.n
        obj.int_p = obj.int_p + alpha(i)*obj.int_pset(i,:);
      end
      obj.int_p = obj.int_p - b*obj.int_pset(obj.n+1,:);

      % get the roots
      r = roots(obj.int_p);
      
      % get logical index of real roots
      rr_lx = abs(imag(r)) <= 1e-10;
      
      % get the largest real root
      rr_max = max(real(r(rr_lx)));
      
      % the final s
      if isempty(rr_max) || rr_max <= 0.0
        s_int = [];
      else
        s_int = 1 / rr_max;
      end
      
    end
    
    function s_bnd = bound(obj,c)
      %bound  compute step sizes to satisfy trust region constraint
      %
      % Given a bound, this method computes step size values such that:
      %
      %   w(s)'*w(s) == c.
      %
      % This equation is always solveable, because ||w(s)|| is an increasing
      % function in s.
      %
      % Input:
      %   c = vector of constraint values
      %
      % Output:
      %   s_bnd = s values corresponding to c
      %
      
      % prepare polynomial data if needed
      % this only need be done once for a given H
      if ~obj.bnd_flg
        obj.bnd_pset = zeros(obj.n+1,2*obj.n+1);
        lx = true(1,obj.n);
        for i = 1:obj.n
          lx(i) = false;
          obj.bnd_pset(i,3:2*obj.n+1) = poly(-[obj.v_hat(lx); obj.v_hat(lx)]);
          lx(i) = true;
        end
        obj.bnd_pset(obj.n+1,:) = poly(-[obj.v_hat; obj.v_hat]);
        obj.bnd_flg = true;
        obj.bnd_p = zeros(1,2*obj.n+1);
      end
      
      % orient and get length of input
      c = c(:)';
      m = length(c);
      s_bnd = zeros(1,m);
      
      for i = 1:m
        
        % handle the non-positive case
        if c(i) <= 0.0
          s_bnd(i) = 0;
          continue;
        end
        
        % compute coefficients
        alpha = obj.utg .* obj.utg;
        
        % construct polynomial for intersection
        obj.bnd_p(:) = 0;
        for j = 1:obj.n
          obj.bnd_p = obj.bnd_p + alpha(j)*obj.bnd_pset(j,:);
        end
        obj.bnd_p = obj.bnd_p - c(i)*obj.bnd_pset(obj.n+1,:);

        % get the roots
        r = roots(obj.bnd_p);
      
        % get logical index of real roots
        rr_lx = abs(imag(r)) <= 1e-10;
      
        % get the largest real root
        rr_max = max(real(r(rr_lx)));
        
        s_bnd(i) = 1/rr_max;
        
      end
      
    end
    
  end
  
end
