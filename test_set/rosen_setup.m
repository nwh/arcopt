function [x0 lb ub n func] = rosen_setup(n)

n = 2;
x0 = [1.2 -1]';
lb = [-5 -5]';
ub = [.5 5]';
func = @rosen;