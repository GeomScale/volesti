function [conv, ratio, too_few] = check_convergence_for_ball(nu, X, r0, lb, ub, lastball)

    conv = false;
    too_few = false;
    
    N = size(X, 2);
    
    mm = N / nu;
    countsIn = 0;
    ratios = [];
    
    for i=1:N
        
        x = X(:,i);
        
        if (norm(x) <= r0)
            countsIn = countsIn + 1;
        end
        
        if (mod(i, mm) == 0)
            ratios = [ratios countsIn/mm];
            countsIn = 0;
        end
    end
    
    ratio = mean(ratios);
    t_a = tinv(0.90, nu);
    rs = std(ratios);
    T = rs * (t_a / sqrt(nu));
    
    if (ratio > lb + T)
        if (lastball)
            conv = true;
            return
        end
        
        if (ratio < ub + T)
            conv = true;
            return
        end
        return
    end
    too_few = true;

end

