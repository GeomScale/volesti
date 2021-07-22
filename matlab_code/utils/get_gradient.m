function [s] = get_gradient(q)
    
    s = nchoosek(q, 2);
    s = s(:,1) .* s(:,2);
    s = s / norm(s);
    
end