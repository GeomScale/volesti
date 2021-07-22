function [Rs, ratios] = get_sequence_of_balls(m, J, L, radius, N, nu, lb, ub)
    
    Rs = [];
    ratios = [];
    n = sum(sum(J>0));

    [r0, ratio0] = get_smallest_ball(m, J, radius, 1200, nu, lb, ub); %ready
    %r0
    %ratio0
    
    %n = sum(sum(J>0));
    Ntot = N * nu;
    
    rmax = Inf;
    rad = Inf;
    while (true)
        
        if (isinf(rmax))
            [points] = billiard_walk_intersection(m, J, L, 4*sqrt(n), Ntot);
        else
            L2 = min(L, rmax);
            [points] = billiard_walk_intersection(m, J, L2, rmax, Ntot);
        end
        
        [conv, ratio] = check_convergence_for_ball(nu, points, r0, lb, ub, true); %ready
        
        if (conv)
            ratios = [ratios ratio];
            Rs = [Rs r0];
            ratios = [ratios ratio0];
            return
        end
        
        if (isinf(rmax))
            rmax = max(sqrt(sum(points.^2,1)));
        else
            rmax = rad;
        end
        
        [rad, ratio] = get_next_ball(points, nu, rmax, r0, lb, ub); %ready
        Rs = [Rs rad];
        ratios = [ratios ratio];
        rmax = rad;
    end
    
    %Rs = [Rs r0];
    %ratios = [ratios ratio0];
    
end