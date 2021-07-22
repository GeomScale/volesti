function [conv, ratio, too_few] = check_convergence_in_spectra(m, J, X, nu, lb, ub)

    too_few = false;
    conv = false;
    
    %n = sum(sum(J>0));
    
    [upper, ~] = initialize_sampler_vol(m);
    A = eye(m);
    
    N = size(X, 2);
    mm = N / nu;
    countsIn = 0;
    ratios = [];
    
    for i=1:N
        
        A(J) = X(:,i);
        q = triu(A',1);
        A(upper) = q(upper);
        
        if (min(eig(A)) > 0)
            countsIn = countsIn + 1;
        end
        
        if (mod(i, mm) == 0)
            ratios = [ratios countsIn/mm];
            countsIn = 0;
            if (length(ratios) > 1)
                
                ratio = mean(ratios);
                t_a = tinv(0.995, nu);
                rs = std(ratios);
                T = rs * (t_a / sqrt(length(ratios)));
                
            
                %[h, ~] = ttest(ratios, lb, 'Alpha', 0.005, 'Tail', 'right');
                if (ratio + T < lb)
                    too_few = true;
                    conv = false;
                    return
                end
            
                %[h, ~] = ttest(ratios, ub, 'Alpha', 0.005, 'Tail', 'left');
                if (ratio - T > ub)
                    conv = false;
                    return
                end
            end
        end
    end
    
    ratio = mean(ratios);
    t_a = tinv(0.95,nu);
    rs = std(ratios);
    T = rs * (t_a / sqrt(nu));
    
    if (ratio > lb + T)
        if (ratio < ub - T)
            conv = true;
            return
        end
        conv = false;
        return
    end
    too_few = true;
    conv = false;

end

