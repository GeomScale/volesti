function [log_vol, Rs, ratios_out, total_points] = volume(m, e, W, J)
    
    n = sum(sum(J>0));
    nu = 10;
    N = 300;
    N_nu = nu * N;
    lb = 0.10;
    ub = 0.15;
    ratios_out = [];

    radius = get_inscribed_radius(m, J);

    L = get_billiard_L(m, 1000)/2;
    
    [Rs, ratios] = get_sequence_of_balls(m, J, L, radius, N, nu, lb, ub);
    
    log_vol = log(pi) * (n/2) + n * log(Rs(end)) - log_gamma_function(n / 2 + 1);
    
    mm = length(Rs) + 1;
    er0 = e / (2.0 * sqrt(mm));
    er1 = (e * sqrt(4.0 * mm - 1)) / (2.0 * sqrt(mm));
    
    ratio = estimate_ratio_first(m, J, Rs(end), ratios(end), er0, W, 1200); 
    ratio_small = ratio;
    log_vol = log_vol + log(ratio);
    
    er1 = er1 / sqrt(mm - 1.0);
    total_points = 0;
    
    if (ratios(1) < 0.9999)
        [ratio, points, numpoints] = estimate_ratio_spectra_ball(m, J, L, Rs(1), [], ratios(1), er1, W, N_nu);
        total_points = total_points + numpoints;
        ratios_out = [ratios_out ratio];
        log_vol = log_vol - log(ratio);
    end
    
    
    for i=1:(mm-2)
        
        L2 = min(Rs(i), L);
        [ratio, points, numpoints] = estimate_ratio(m, J, L2, Rs(i), Rs(i+1), points, ratios(i+1), er1, W, N_nu); %ready
        %ratio
        total_points = total_points + numpoints;
        ratios_out = [ratios_out ratio];
        log_vol = log_vol - log(ratio);
        
    end
    ratios_out = [ratios_out ratio_small];

end



