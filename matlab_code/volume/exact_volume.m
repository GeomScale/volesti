function log_vol = exact_volume(m)
    
    a = (m+1)/2;
    log_vol = -m * log_gamma_function((m+1)/2);

    log_vol = log_vol + ((m*(m-1))/4)*log(pi);
    
    for i=1:m
        log_vol = log_vol + log_gamma_function(a - ((i-1)/2));
    end
end

