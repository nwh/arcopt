function [x0 lb ub n func] = genrose_setup(n)

x0 = ones(n,1);
x0(1) = -1.2;
if n >= 3
  x0(3) = -1.2;
end

x = ones(n,1);

lb = -100*ones(n,1);
ub = 100*ones(n,1);

for i = 1:2:n
  lb(i) = x(i) + 0.1;
  ub(i) = x(i) + 1.1;
end

func = @genrose;