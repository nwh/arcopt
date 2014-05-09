function [f g H] = rosen(x)
% rosen rosenbrock test function
%
% 2010-02-17 (nwh) created rosen.m

f = (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;

if ( nargout > 1 )
  % return the gradient
  g = [0;0];
  g(1) = 400*x(1)^3-2*x(1)*(200*x(2)-1)-2;
  g(2) = 200*(x(2)-x(1)^2);
  
  if ( nargout > 2 )
    % return the Hessian
    H = zeros(2,2);
    H(1,1) = 1200*x(1)^2 - 2*(200*x(2)-1);
    H(1,2) = -400*x(1);
    H(2,1) = -400*x(1);
    H(2,2) = 200;
  end
end
end