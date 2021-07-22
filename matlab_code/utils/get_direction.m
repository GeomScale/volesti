function [v] = get_direction(n)

    v = randn(n,1);
    v = v / norm(v);

end