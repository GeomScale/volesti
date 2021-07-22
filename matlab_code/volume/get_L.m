function L = get_L(m, J, radius)

    n = m^2 / 2 - m/2;
    L = 5 * sqrt(n) * radius;
    
    points = billiard_walk_intersection(m, J, L, 10^10, 10 * (ceil(log(n))*10));
    
    L = max([L, 2 * max(sqrt(sum(points.^2, 1)))]);

end