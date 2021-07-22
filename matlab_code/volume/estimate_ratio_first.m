function ratio = estimate_ratio_first(m, J, R, ratio, error, W, N_nu)

    n = sum(sum(J>0));
    
    [upper, ~] = initialize_sampler(m);
    A = eye(m);

    conv = false;
    
    ratio_parameters = initialize_parameters(W, N_nu, ratio);
    
    while (~conv)
        
        x = randsphere(1, n, R)';
        A(J) = x;
        q = triu(A',1);
        A(upper) = q(upper);
        
        is_in = false;
        if (min(eig(A)) > 0)
            is_in = true;
        end
        
        [conv, ratio, ratio_parameters] = update_window(ratio_parameters, is_in, error);
       
    end
    
    %ratio = ratio_parameters.count_in / ratio_parameters.tot_count;
    
end