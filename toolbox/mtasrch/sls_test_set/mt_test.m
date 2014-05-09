
f = {@(x)mt1(x); @(x)mt2(x); @(x)mt3(x); @(x)mt4(x)};
b = {[0 16]; [0 2]; [0 2]; [0 1]};

n = 200;
y = zeros(n,2);

for i = 1:length(f)
  
  x = linspace(b{i}(1),b{i}(2),n);
  
  for j = 1:n
    [y(j,1) y(j,2)]= f{i}(x(j));
  end
  
  figure(1)
  subplot(2,1,1)
  plot(x,y(:,1))
  subplot(2,1,2)
  plot(x,y(:,2))
  pause
end