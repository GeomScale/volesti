addpath('./utils')
addpath('./volume')

%size of correlation matrices
m=8;
%initialize billiard walk sampler
[upper, lower] = initialize_sampler(m);
J = lower;
%estimate the logarithm of the volume
[log_vol, Rs, ratios_out, tot_points] = volume(m, 0.1, 1500, J);

esti_vol = exp(log_vol);
disp(strcat('estimated volume = ', num2str(esti_vol)));

exact_vol = exp(exact_volume(m));
disp(strcat('exact volume = ', num2str(exact_vol)));

rel_error = abs(exp(log_vol) - exp(exact_volume(m)))/ exp(exact_volume(m));
disp(strcat('relative error = ', num2str(rel_error)));
