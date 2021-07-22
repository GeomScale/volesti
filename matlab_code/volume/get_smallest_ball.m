function [r, ratio] = get_smallest_ball(m, J, radius, N, nu, lb, ub)

    n = sum(sum(J>0));
    tolerance = 0.00000000001;
    
    sqrt_n = sqrt(n);
    
    rad1 = radius;
    rmax = 2 * sqrt_n * radius;
    conv = false;
    
    while (~conv)
      
        X = randsphere(N, n, rmax)';
        %size(X)
        
        [conv, ratio, too_few] = check_convergence_in_spectra(m, J, X, nu, lb, ub);

        if (conv)
            r = rmax;
            return
        end
        if (too_few)
            break
        end
        rad1 = rmax;
        rmax = rmax + 2 * sqrt_n * radius;
    end
    
    rad0 = rad1;
    rad_m = rmax;
    while (true)
        
        rad_med = 0.5 * (rad1 + rmax);
        
        X = randsphere(N, n, rad_med)';
        [conv, ratio, too_few] = check_convergence_in_spectra(m, J, X, nu, lb, ub);
        
        if (conv)
            r = rad_med;
            return
        end
        
        if (too_few)
            rmax = rad_med;
        else
            rad1 = rad_med;
        end
        
        if (rmax - rad1 < tolerance)
            rad1 = rad0;
            rmax = rad_m;
        end
        
    end

end