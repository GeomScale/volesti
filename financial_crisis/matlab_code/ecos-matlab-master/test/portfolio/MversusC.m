% compares values written by Matlab and C-version of ECOS

% run Matlab version
conelp_solver;

% run C version
ecos_solver;

ok = 1;

threshold = 1e-14;

fprintf('* ECOS MATLAB VERSUS C CHECKER (threshold: %4.2e) *\n', threshold);
fprintf('-------------------------------------------------------\n');

%% Test permutation vector
ok = ok && compareMatrices('M_P.txt', 'P.txt', threshold);
fprintf('-------------------------------------------------------\n');

%% Test initialization
ok = ok && compareMatrices('M_x_init.txt', 'x_init.txt', threshold);
ok = ok && compareMatrices('M_y_init.txt', 'y_init.txt', threshold);
ok = ok && compareMatrices('M_z_init.txt', 'z_init.txt', threshold);
ok = ok && compareMatrices('M_s_init.txt', 's_init.txt', threshold);
fprintf('-------------------------------------------------------\n');