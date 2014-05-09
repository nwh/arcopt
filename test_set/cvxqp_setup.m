function [x0 lb ub n func] = cvxqp_setup(n);

n = 20;
x0 = zeros(n,1);
lb = -2*ones(n,1);
ub = 2*ones(n,1);
func = @cvxqp;