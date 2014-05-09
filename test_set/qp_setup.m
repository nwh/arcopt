function [x0 lb ub n func] = qp_setup(n)

n = 10;
x0 = zeros(n,1);
lb = -10*ones(n,1);
ub = 10*ones(n,1);
func = @qp;