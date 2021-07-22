function L = get_billiard_L(m, N)

    n = m^2 / 2 - m/2;
    
    [upper, lower] = initialize_sampler(m);
    x = zeros(n, 1);
    A = eye(m);
    B = zeros(m);
    
    bc = ones(2 * n, 1);
    L = 0;
    
    for j = 1:N
            v = get_direction(n);
            
            %A(lower) = x;
            %q = triu(A',1);
            %A(upper) = q(upper);
                
            B(lower) = v;
            q = triu(B',1);
            B(upper) = q(upper);
            
            % compute the intersection of the line x + l*v with the
            % boundary of the convex set that contains the PSD matrices
            % with ones in the diagonal
            eigenvalues = eig(B, -A);
            l_max = 1 / max(eigenvalues);
            l_min = 1 / min(eigenvalues);
            
            % compute the intersection of the line x + l*v with the
            % boundary of the hypercube [-1, 1]^n
            lambdas = [v; -v] ./ (bc - [x; -x]);
            l_max_temp = 1 / max(lambdas);
            l_min_temp = 1 / min(lambdas);
                
            % decide which boundary the line hits first
            if (l_max_temp < l_max)
                l_max = l_max_temp;
            end
            if (l_min_temp > l_min)
                l_min = l_min_temp;
            end
            
            if (L < (l_max + l_min))
                L = l_max + l_min;
            end
            
            % pick a uniformly distributed point from the segment
            %lambda = l_min + rand * (l_max - l_min);
            
            % update the current point of the random walk
            %x = x + lambda * v;
            
    end
    
    L = 4 * L;

end
