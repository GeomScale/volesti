function y = log_gamma_function(x) 

    if (x <= 100)
        y = log(gamma(x));
        return
    end
    
    y = log(x - 1) + log_gamma_function(x - 1);
end