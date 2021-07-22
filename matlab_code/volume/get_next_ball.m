function [rad, ratio] = get_next_ball(X, nu, rmax, r0, lb, ub)

    rad0 = r0;
    rad_m = rmax;
    tolerance = 0.00000000001;
   
    while (true)
        
        rad_med = 0.5 * (r0 + rmax);
        
        [conv, ratio, too_few] = check_convergence_for_ball(nu, X, rad_med, lb, ub, false);
    
        if (conv)
            rad = rad_med;
            return
        end
        
        if (too_few)
            r0 = rad_med;
        else
            rmax = rad_med;
        end
        
        if (rmax - r0 < tolerance)
            r0 = rad0;
            rmax = rad_m;
        end
        
    end

end


