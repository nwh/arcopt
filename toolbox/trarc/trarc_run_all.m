%trarc_run_all  run all trarc scripts

% run the basic test
trarc_run;

% run the intersection test
trarc_run_int;

% run the bound test
trarc_run_bnd;

% plot box intersection
% use figure(1)
trarc_plt_box(1);

% plot trust region intersection
% use figure(2)
trarc_plt_cir(2);

% run the exp tests
trarc_run_exp;
