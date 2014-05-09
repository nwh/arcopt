function [f g] = var28(x,lambda)

if nargin < 2
  lambda = -3.4;
end

% orient vector
x = x(:);

% get parameters
n = length(x);
h = 1/(n+1);

% pad vector
y = [0; x; 0];

% compute differences and exponentials
dy = diff(y);
expy = exp(y);
dexpy = diff(expy);

% compute function value
f = 2/h*sum(-x.*dy(2:end))+2*lambda*h*sum(dexpy./dy);

% compute gradient
g1 = diff(y,2);
g2 = ( exp(x) .* dy(1:end-1) - dexpy(1:end-1)) ./ (dy(1:end-1).^2);
g3 = (-exp(x) .* dy(2:end)   + dexpy(2:end)  ) ./ (dy(2:end)  .^2);

g = -2/h*g1+2*lambda*h*(g2+g3);
