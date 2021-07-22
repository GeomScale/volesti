function [conv, val, ratio_parameters] = update_window(ratio_parameters, is_in, error)

    conv = false;
    
    if (is_in)
        ratio_parameters.count_in = ratio_parameters.count_in + 1;
    end
    
    ratio_parameters.tot_count = ratio_parameters.tot_count + 1;
    val = ratio_parameters.count_in / ratio_parameters.tot_count;
    ratio_parameters.last_W(ratio_parameters.index) = val;

    if (val <= ratio_parameters.min_val)
        ratio_parameters.min_val = val;
        ratio_parameters.min_index = ratio_parameters.index;
    elseif (ratio_parameters.min_index == ratio_parameters.index)
        [ratio_parameters.min_val, ratio_parameters.min_index] = min(ratio_parameters.last_W);
        %ratio_parameters.minmaxIt = std::min_element(ratio_parameters.last_W.begin(), ratio_parameters.last_W.end());
        %ratio_parameters.min_val = (*ratio_parameters.minmaxIt);
        %ratio_parameters.min_index = std::distance(ratio_parameters.last_W.begin(), ratio_parameters.minmaxIt);
    end

    if (val >= ratio_parameters.max_val)
        ratio_parameters.max_val = val;
        ratio_parameters.max_index = ratio_parameters.index;
    elseif (ratio_parameters.max_index == ratio_parameters.index)
        [ratio_parameters.max_val, ratio_parameters.max_index] = max(ratio_parameters.last_W);
        %ratio_parameters.minmaxIt = std::max_element(ratio_parameters.last_W.begin(), ratio_parameters.last_W.end());
        %ratio_parameters.max_val = (*ratio_parameters.minmaxIt);
        %ratio_parameters.max_index = std::distance(ratio_parameters.last_W.begin(), ratio_parameters.minmaxIt);
    end
    
    %ratio_parameters.last_W
    %ratio_parameters.min_val
    %ratio_parameters.max_val
    
    %ratio_parameters.last_W
    %ratio_parameters.max_val
    %ratio_parameters.max_index
    
    %ratio_parameters.min_val
    %ratio_parameters.min_index
    
    %error
    if ( (ratio_parameters.max_val - ratio_parameters.min_val) / ratio_parameters.max_val <= error/2 )
        %(ratio_parameters.max_val - ratio_parameters.min_val) / ratio_parameters.max_val
        %val
        %ratio_parameters.tot_count
        conv = true;
        return
    end
       

    ratio_parameters.index = mod(ratio_parameters.index, ratio_parameters.W) + 1;
    %if (ratio_parameters.index == ratio_parameters.W + 1) 
    %    ratio_parameters.index = 1;
    %end
    

end