function [phi dphi] = mt2(alpha,param)
% mt2: More and Thuente test function 2
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
%
% Usage:
%   [phi dphi] = mt2(a);
%
%   param.beta = 1;
%   [phi dphi] = mt2(a,param)
%
% 2010-02-06 (nwh) created test function

if nargin < 2
  beta = 0.004;
else
  beta = param.beta;
end

phi = (alpha+beta)^5 - 2*(alpha+beta)^4;
dphi = 5*(alpha+beta)^4 - 8*(alpha+beta)^3;
