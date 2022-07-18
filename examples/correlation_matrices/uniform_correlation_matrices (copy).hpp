template<typename WalkTypePolicy,
    typename PointType,
    typename PointList,
    typename RNGType
    >
void uniform_correlation_sampling(CorreSpectra<PointType> &P,
                                    PointList &randPoints,
                                    RNGType &rng,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    const PointType &starting_point,
                                    unsigned int const& nburns){
    typedef typename WalkTypePolicy::template Walk <CorreSpectra<PointType>, RNGType> walk;
    PushBackWalkPolicy push_back_policy;
    typedef RandomPointGenerator<walk> RandomPointGenerator;
    
    PointType p = starting_point;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walkL, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, num_points, walkL, randPoints,
                                push_back_policy, rng);
}
f = @(x, y) exp(-0.3*x'*x);
logf_neg = @(x, y) 0.3*x'*x;
samples = ReHMC(m, N, walk_length, 0.03, 1, logf_neg);
function [samples] = ReHMC(m, N, W, step, L, f, f_utils, T)
% Reflective Hamiltonian Monte Carlo to sample from f truncated to the set
% of the mxm correlation matrices
    
    n = m^2 / 2 - m/2;
    w = ceil(L / step);
    [upper, lower] = initialize_sampler(m);
    x0 = zeros(n, 1);
    A = eye(m);
    B = zeros(m);
    
    bc = ones(2 * n, 1);
    pos = 1:(2*n);
    pos = mod(pos, n);
    pos(pos==0) = n;
    samples = cell(N, 1);
    
    f_utils.T = T;
    
    h = waitbar(0,'Computing samples...');
    for i = 1:N
        for j = 1:W
            T = step;
            v = randn(n, 1);
            x = x0;
            grad_x = esti_grad(f, f_utils, x);
            h1 = f(x, f_utils) + 0.5 * (v' * v);
            for k = 1:w
                
                v = v - (step/2) * grad_x;
                
                while(true)
                
                    [Q, eigenvalues] = eig(B, -A);
                    [max_eig, pos_max_eig] = max(diag(eigenvalues));
                    l_max = 1 / max_eig;
                
                    lambdas = [v; -v] ./ (bc - [x; -x]);
                    [l_max_temp, pos_max] = max(lambdas);
                    l_max_temp = 1 / l_max_temp;
                
                    
                    if (l_max_temp < l_max)
                
                       
                        if (T <= l_max)
                            x = x + T * v;
                            break;
                        end
                        
                        lambda = 0.995 * l_max;
                        s = get_gradient(Q(:, pos_max_eig));
                        % update the current point of the random walk
                        x = x + lambda * v;
                        %reflect the ray
                        v = v - (2*(v'*s))*s;
                    
                    T = T - lambda;
                end
                grad_x = esti_grad(f, f_utils, x);
                v = v - (step/2) * grad_x;
            end
            h2 = f(x, f_utils) + 0.5 * (v' * v);
            u = rand;
            if (u < exp(h1-h2))
                x0 = x;
            end        
        end
        A(lower) = x0;
        q = triu(A',1);
        A(upper) = q(upper); 
        
        samples{i} = A;
        waitbar(i/N);
    end    
    close(h);
end

