function [ratio, points, numpoints] = estimate_ratio_spectra_ball(m, J, L, R, p, ratio, error, W, N_nu)

    n = sum(sum(J>0));
    points = p;
    
    [upper, ~] = initialize_sampler(m);
    %x = zeros(n, 1);
    A = eye(m);
    B = zeros(m);
    
    bc = ones(2 * n, 1);
    pos = 1:(2*n);
    pos = mod(pos, n);
    pos(pos==0) = n;
    
    conv = false;
    ratio_parameters = initialize_parameters(W, N_nu, ratio);
    
    while (~conv)
        
        for jj=1:3
        T = -log(rand) * L;
        v = get_direction(n);
            
        rho = 0;
        
        while (true)
            
            %A(lower) = x;
            %q = triu(A',1);
            %A(upper) = q(upper);
                
            B(J) = v;
            q = triu(B',1);
            B(upper) = q(upper);
            
            % compute the intersection of the line x + l*v with the
            % boundary of the convex set that contains the PSD matrices
            % with ones in the diagonal
            [Q, eigenvalues] = eig(B, -A);
            [max_eig, pos_max_eig] = max(diag(eigenvalues));
            l_max = 1 / max_eig;
            
            x = A(J);
            
            % compute the intersection of the line x + l*v with the
            % boundary of the hypercube [-1, 1]^n
            lambdas = [v; -v] ./ (bc - [x; -x]);
            [l_max_temp, pos_max] = max(lambdas);
            l_max_temp = 1 / l_max_temp;
            
            
            
            [l_max, lmax_ind] = min([l_max l_max_temp]);
                
            % decide which boundary the line hits first
            if (lmax_ind == 2)
                %l_max = l_max_temp;
                % pick a uniformly distributed point from the segment
                lambda = 0.995 * l_max;
            
                % update the current point of the random walk
                if (T <= l_max)
                    %x = x + T * v;
                    A = A + T * B;
                    break;
                end
                A = A + lambda * B;
                %x = x + lambda * v;
                    
                %reflevt the ray
                %v = v - (2*AA(pos_max, :)*v)*AA(pos_max, :)'
                p = pos(pos_max);
                v(p) = -v(p);
            elseif (lmax_ind == 1)
                % pick a uniformly distributed point from the segment
                lambda = 0.995 * l_max;
                if (T <= l_max)
                    %x = x + T * v;
                    A = A + T* B;
                    break;
                end
                s = get_gradient(Q(:, pos_max_eig));
                % update the current point of the random walk
                %x = x + lambda * v;
                A = A + lambda * B;
                %reflect the ray
                v = v - (2*(v'*s))*s;
            end
            rho = rho + 1;
            T = T - lambda;
        end
        end
        %A(lower) = x;
        %q = triu(A',1);
        %A(upper) = q(upper);
        x = A(J);
        
        is_in = false;
        if (norm(x) <= R)
            is_in = true;
        end
        
        [conv, ratio, ratio_parameters] = update_window(ratio_parameters, is_in, error);
        
    end
    numpoints = ratio_parameters.tot_count;
    %ratio = ratio_parameters.count_in / ratio_parameters.tot_count;
    

end