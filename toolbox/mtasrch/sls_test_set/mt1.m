function [phi dphi] = mt1(alpha,param)
% mt1: More and Thuente test function 1
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
%   [phi dphi] = mt1(a);
%
%   param.beta = 1;
%   [phi dphi] = mt1(a,param)
%
% 2010-02-06 (nwh) created test function

if nargin < 2
  beta = 2;
else
  beta = param.beta;
end

phi = -alpha/(alpha^2+beta);
dphi = (alpha^2-beta)/(alpha^2+beta)^2;
