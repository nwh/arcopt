% setup script for arcopt

function load_arcopt

home_dir = getenv('HOME');
arcopt_dir = [home_dir '/Dropbox/code/arcopt'];
addpath(arcopt_dir);
addpath([arcopt_dir '/toolbox/mtasrch']);
addpath([arcopt_dir '/toolbox/trarc']);
addpath([arcopt_dir '/toolbox/lusol_mex']);
addpath([arcopt_dir '/toolbox/expand']);

end
