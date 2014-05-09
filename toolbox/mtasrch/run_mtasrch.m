% run_mtas.m
%
% Test script for mtas, or Brain Dead More and Thuente step length
% selection.
%
% 2010-02-06 (nwh) created test script
% 2010-02-16 (nwh) modified to use mtas
% 2010-02-25 (nwh) handle mtas output
%

% add the test problems to the path
addpath sls_test_set

% set parameters for test problem mt4
p2.beta1 = 0.01;
p2.beta2 = 0.001;
p3.beta1 = 0.001;
p3.beta2 = 0.01;

% construct test set table
% data fields:
%  name func a0 amin amax ftol gtol itol
test_set = {
  'mt1' @mt1 1e-3 0 2000 0.001 0.1 1e-6
  'mt1' @mt1 1e-1 0 2000 0.001 0.1 1e-6
  'mt1' @mt1 1    0 2000 0.001 0.1 1e-6
  'mt1' @mt1 1e1  0 2000 0.001 0.1 1e-6
  'mt1' @mt1 1e3  0 2000 0.001 0.1 1e-6
  'mt2' @mt2 1e-3 0 2000 0.1 0.1 1e-10
  'mt2' @mt2 1e-1 0 2000 0.1 0.1 1e-10
  'mt2' @mt2 1    0 2000 0.1 0.1 1e-10
  'mt2' @mt2 1e1  0 2000 0.1 0.1 1e-10
  'mt2' @mt2 1e3  0 2000 0.1 0.1 1e-10
  'mt3' @mt3 1e-3 0 2000 0.1 0.1 1e-10
  'mt3' @mt3 1e-1 0 2000 0.1 0.1 1e-10
  'mt3' @mt3 1    0 2000 0.1 0.1 1e-10
  'mt3' @mt3 1e1  0 2000 0.1 0.1 1e-10
  'mt3' @mt3 1e3  0 2000 0.1 0.1 1e-10
  'mt4' @mt4 1e-3 0 2000 0.001 0.001 1e-10
  'mt4' @mt4 1e-1 0 2000 0.001 0.001 1e-10
  'mt4' @mt4 1    0 2000 0.001 0.001 1e-10
  'mt4' @mt4 1e1  0 2000 0.001 0.001 1e-10
  'mt4' @mt4 1e3  0 2000 0.001 0.001 1e-10
  'mt5' @(x)mt4(x,p2) 1e-3 0 2000 0.001 0.001 1e-10
  'mt5' @(x)mt4(x,p2) 1e-1 0 2000 0.001 0.001 1e-10
  'mt5' @(x)mt4(x,p2) 1    0 2000 0.001 0.001 1e-10
  'mt5' @(x)mt4(x,p2) 1e1  0 2000 0.001 0.001 1e-10
  'mt5' @(x)mt4(x,p2) 1e3  0 2000 0.001 0.001 1e-10
  'mt6' @(x)mt4(x,p3) 1e-3 0 2000 0.001 0.001 1e-10
  'mt6' @(x)mt4(x,p3) 1e-1 0 2000 0.001 0.001 1e-10
  'mt6' @(x)mt4(x,p3) 1    0 2000 0.001 0.001 1e-10
  'mt6' @(x)mt4(x,p3) 1e1  0 2000 0.001 0.001 1e-10
  'mt6' @(x)mt4(x,p3) 1e3  0 2000 0.001 0.001 1e-10
};

[num_tests junk] = size(test_set);
tests = 1:num_tests;

% prepare string formats for outputting
header_str = '%2s %8s %8s %2s %3s %8s %8s\n';
output_str = '%2d %8s %8.3g %2d %3d %8.3g %8.2g\n';
fprintf(header_str,'i','prob','a0','tf','nf','a*','dphi(a*)')
fprintf(char(45*ones(45,1)))
fprintf('\n')

% carry out all tests
for i = tests
  name = test_set{i,1};
  fcn = test_set{i,2};
  a0 = test_set{i,3};
  stpmin = test_set{i,4};
  stpmax = test_set{i,5};
  ftol = test_set{i,6};
  gtol = test_set{i,7};
  xtol = test_set{i,8};

  s = 1;
  n = 1;
  x = 0;
  [f g] = fcn(x);
  maxfev = 100;
  
  % open output file
  fid = fopen(['output/' name '_' num2str(i) '.txt'],'w');
  
  % construct anonymous function for arc
  arc = @(stp) line_arc(stp,s);
  
  % call bdmt
  [x f g stp info nfev] = mtasrch(fcn,x,f,g,arc,a0, ...
                                  0,ftol,gtol,xtol, ...
                                  stpmin,stpmax,maxfev,fid);
  
  fclose(fid);
  fprintf(output_str,i,name,a0,info,nfev,stp,g);
end
