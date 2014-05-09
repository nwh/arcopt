% setup script for arcopt

function arcopt_setup

load_mcute
load_xunit
load_arcopt

home_dir = getenv('HOME');
arcopt_dir = [home_dir '/Dropbox/code/arcopt'];
addpath(arcopt_dir);
addpath([arcopt_dir '/test_set']);

end
