function parameters = initialize_parameters(W, N, ratio)

    parameters = struct;
    
    parameters.min_val = -Inf;
    parameters.max_val = Inf;
    parameters.min_index = W;
    parameters.max_index = W;
    parameters.W = W;
    parameters.index =1;
    parameters.tot_count = N;
    parameters.count_in = floor(N * ratio);
    %parameters.tot_count = 0;
    %parameters.count_in = 0;
    parameters.last_W = zeros(1, W);

end