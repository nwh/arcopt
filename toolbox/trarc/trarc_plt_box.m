%trarc_plt_box  plot a trust region arc in a box

function trarc_plt_box(fig)
  
  % set optional value for figure
  if nargin < 1 || isempty(fig)
    fig = 1;
  end
  
  % set the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % generate data
  n = 2;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = randn(n,1);
  
  % set box size
  c = 2;
  
  % get arc
  arc = trarc(H,g);
  
  % construct box constraint
  m = 2*n;
  A = [eye(n); -eye(n)];
  b = -c*ones(m,1);
  
  % compute the intersection
  s_max = 100;
  s_ix = 0;
  for i = 1:m
    s_int = arc.intersect(A(i,:),b(i));
    if ~isempty(s_int) && s_int < s_max
      s_max = s_int;
      s_ix = i;
    end
  end
  
  % plot
  figure(fig); clf; hold on
  cp = c + 0.1*c;
  s = linspace(0,s_max,100);
  w = arc.sol_s(s);
  plot(w(1,:),w(2,:));
  axis([-cp cp -cp cp])
  axis square
  
  % plot a box
  plot([-c c c -c -c],[-c -c c c -c],'k')
  
  % plot circle at end of arc
  plot(w(1,end),w(2,end),'ro')
  
end
