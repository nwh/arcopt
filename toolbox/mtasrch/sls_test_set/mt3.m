function [phi dphi] = mt3(alpha,param)
% mt3: More and Thuente test function 3
%
% Input:
%   alpha = step length
%   param = optional parameter struct
%
% Output:
%   phi = function value
%   dphi = first derivative
%
% param fields:
%   beta
%   l
%
% Usage:
%   [phi dphi] = mt3(a);
%
%   param.beta = 1;
%   param.l = 3;
%   [phi dphi] = mt3(a,param)
%
% 2010-02-06 (nwh) created test function

if nargin < 2
  beta = 0.01;
  l = 39;
else
  beta = param.beta;
  l = param.l;
end

if alpha <= 1-beta
  phi0 = 1-alpha;
  dphi0 = -1;
elseif alpha >= 1+beta
  phi0 = alpha-1;
  dphi0 = 1;
else
  phi0 = (alpha-1)^2/2/beta + beta/2;
  dphi0 = (alpha-1)/beta;
end

phi = phi0 + 2*(1-beta)/l/pi*sin(alpha*pi*l/2);
dphi = dphi0 + (1-beta)*cos(alpha*pi*l/2);
