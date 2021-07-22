function [upper, lower] = initialize_sampler(m)
    
    At = zeros(m, m);
    upper  = (1:size(At,1)).' < (1:size(At,2));
    lower = (1:size(At,1)).' > (1:size(At,2));
    
end