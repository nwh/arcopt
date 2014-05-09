%trarc_plt_cir  plot a 2d test of trarc.bound

function trarc_plt_cir(fig)
  
  % set optional value for figure
  if nargin < 1 || isempty(fig)
    fig = 1;
  end
  
  % set the rng
  % RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
  
  % data
  n = 2;
  H = randn(n,n);
  H = 0.5*(H+H');
  g = -randn(n,1);
  
  % tr bound
  c = [1 2 3];
  m = length(c);
  
  % get arc
  arc = trarc(H,g);
  
  % get the step sizes
  s_bnd = arc.bound(c);
  w_bnd = arc.sol_s(s_bnd);
  
  % plot
  
  figure(fig); clf; hold on
  
  s = linspace(0,max(s_bnd),100);
  w = arc.sol_s(s);
  plot(w(1,:),w(2,:));
  c1 = sqrt(max(c))+.1;
  axis([-c1 c1 -c1 c1])
  axis square
  
  % plot circles
  for i = 1:m
    circle([0 0],sqrt(c(i)),300,'k-');
  end
  
  % plot intersections
  plot(w_bnd(1,:),w_bnd(2,:),'ro');
  
end

function H=circle(center,radius,NOP,style)
  %---------------------------------------------------------------------------------------------
  % H=CIRCLE(CENTER,RADIUS,NOP,STYLE)
  % This routine draws a circle with center defined as
  % a vector CENTER, radius as a scaler RADIS. NOP is
  % the number of points on the circle. As to STYLE,
  % use it the same way as you use the rountine PLOT.
  % Since the handle of the object is returned, you
  % use routine SET to get the best result.
  %
  %   Usage Examples,
  %
  %   circle([1,3],3,1000,':');
  %   circle([2,4],2,1000,'--');
  %
  %   Zhenhai Wang <zhenhai@ieee.org>
  %   Version 1.00
  %   December, 2002
  %---------------------------------------------------------------------------------------------
  
  if (nargin <3),
    error('Please see help for INPUT DATA.');
  elseif (nargin==3)
    style='b-';
  end;
  THETA=linspace(0,2*pi,NOP);
  RHO=ones(1,NOP)*radius;
  [X,Y] = pol2cart(THETA,RHO);
  X=X+center(1);
  Y=Y+center(2);
  H=plot(X,Y,style);
  axis square;
  
end