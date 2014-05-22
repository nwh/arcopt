% setup script for arcopt
%
% Usage:
%   >> cd [arcopt directory]
%   >> arcopt_setup('install') % to create arcopt_load file in userpath
%   >> arcopt_setup            % to add arcopt directly to path
%
% After arcopt_setup('install'), the following command will add arcopt to
% Matlab's path:
%   >> arcopt_load
%
function arcopt_setup(install)

if nargin == 1 && strcmp(install,'install')
  % get userpath directory
  mypath = strsplit(userpath,':');
  mypath = mypath{1};
  % create directory if it does not exist
  % TODO: make this better
  [~,~,~] = mkdir(mypath)
  % get arcopt dir
  arcopt_dir = pwd;
  % open file
  fid = fopen([mypath '/load_arcopt.m'],'w');
  % write file
  fprintf(fid,'function load_arcopt\n');
  fprintf(fid,'  addpath(''%s'')\n',arcopt_dir);
  fprintf(fid,'  addpath(''%s'')\n',[arcopt_dir '/toolbox/mtasrch']);
  fprintf(fid,'  addpath(''%s'')\n',[arcopt_dir '/toolbox/trarc']);
  fprintf(fid,'  addpath(''%s'')\n',[arcopt_dir '/toolbox/lusol/matlab']);
  fprintf(fid,'  addpath(''%s'')\n',[arcopt_dir '/toolbox/expand']);
  fprintf(fid,'  addpath(''%s'')\n',[arcopt_dir '/test_set']);
  fprintf(fid,'end\n');
  % close file
  fclose(fid);
  %keyboard
else
  arcopt_dir = pwd;
  addpath(arcopt_dir);
  addpath([arcopt_dir '/toolbox/mtasrch']);
  addpath([arcopt_dir '/toolbox/trarc']);
  addpath([arcopt_dir '/toolbox/lusol/matlab']);
  addpath([arcopt_dir '/toolbox/expand']);
  addpath([arcopt_dir '/test_set']);
end

end
