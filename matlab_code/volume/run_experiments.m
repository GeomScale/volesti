mat_sizes = [8 10 12];
errors = zeros(3, 10);
runtimes = zeros(3, 10);
log_vols = zeros(3, 10);
total_points = zeros(3,10);

for j = 1:length(mat_sizes)
    
    m = mat_sizes(j);
    [upper, lower] = initialize_sampler(m);
    J = lower;

    for i=1:10
        disp([m i])
        tic
        [log_vol, Rs, ratios_out, tot_points] = volume(m, 0.1, 1000*j, J);
        tim = toc;
        log_vols(j, i) = log_vol;
        runtimes(j, i) = tim;
        total_points(j, i) = tot_points;
        errors(j, i) = abs(exp(log_vol) - exp(exact_volume(m)))/ exp(exact_volume(m))
    end
    
end
save('elliptope_volume_results.mat', 'log_vols', 'runtimes', 'errors')

