function radius = get_inscribed_radius(m, J)

    n = sum(sum(J));
    
    [upper, ~] = initialize_sampler_vol(m);
    A = eye(m);
    B = zeros(m);
    
    radius = Inf;
    
    for i=1:n
       
        v = zeros(n,1);
        v(i) = 1;
        
        B(J) = v;
        q = triu(B',1);
        B(upper) = q(upper);
        
        [~, eigenvalues] = eig(B, -A);
        max_eig = max(diag(eigenvalues));
        l_pos = 1 / max_eig;
        
        [~, eigenvalues] = eig(-B, -A);
        max_eig = max(diag(eigenvalues));
        l_neg = 1 / max_eig;
        
        ell = min([l_pos, l_neg]);
        
        if (ell < radius)
            radius = ell;
        end
        
    end
    
    radius = radius / sqrt(n);

end