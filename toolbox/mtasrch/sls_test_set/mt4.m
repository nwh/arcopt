function [phi dphi] = mt4(alpha,param)
% mt4: More and Thuente test function 4
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
%   beta1
%   beta2
%
% Usage:
%   [phi dphi] = mt4(a);
%
%   param.beta1 = 1;
%   param.beta2 = 3;
%   [phi dphi] = mt4(a,param)
%
% 2010-02-06 (nwh) created test function

if nargin < 2
  beta1 = 0.001;
  beta2 = 0.001;
else
  beta1 = param.beta1;
  beta2 = param.beta2;
end

gamma = @(b) sqrt(1+b^2)-b;

phi = gamma(beta1)*sqrt((1-alpha)^2+beta2^2) + ...
      gamma(beta2)*sqrt(alpha^2+beta1^2);
dphi = gamma(beta1)*(alpha-1)/sqrt((1-alpha)^2+beta2^2) + ...
      gamma(beta2)*alpha/sqrt(alpha^2+beta1^2);
