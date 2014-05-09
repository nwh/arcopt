function [x0 lb ub n func] = var28_setup(n)

% parameters
h = 1/(n+1);

% initial point
x0 = 0.1.*(1:n).*h.*(-(1:n));
x0 = x0';

% bounds
lb = -0.2*n*ones(n,1);
ub = 0.2*n*ones(n,1);

func = @var28;