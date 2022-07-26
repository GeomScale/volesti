/// ReHMC walk

template<typename NT, typename Point>
NT Hamiltonian(Point const &p, Point const &v){
    return (p.dot(p) + v.dot(v))/2;
}

template<typename Point>
Point esti_grad(Point const &p){
    return p;
}

template<typename NT, typename ConvexBodyType, typename Point, typename RNGType>
Point ReHMC(ConvexBodyType &P, Point& q, unsigned int const& walkL, RNGType &rng, NT const step){
    int n = P.dimension(), num_leaps = ceil(walkL / step);
    int i, j, it, nreflex = 50*n;
    NT T, h1, h2, lambda, u;
    Point p0 = q, p = q, v;
    Point grad_p;
    
    for(i = 0; i < walkL; ++i){
        T = step; // leap total distance T = rng.sample_urdist() * _Len;
        v = GetDirection<Point>::apply(n, rng, false);
        p = p0;

        h1 = Hamiltonian<NT>(p, v);
        // std::cout << "h1 = " << h1 << std::endl;
        grad_p = esti_grad<Point>(p);

        for(j = 0; j < num_leaps; ++j){            
            v = v - (step/2) * grad_p;

            // it = 0;
            while(true){
            // while(it < nreflex){
                
                auto pbpair = P.line_positive_intersect(p, v);
                if (T <= pbpair.first){
                    p += T * v;
                    break;
                }
                
                lambda = 0.995 * pbpair.first;
                T -= lambda;
                p += lambda * v;
                P.compute_reflection(v, p, pbpair.second);
                // it++;
            }
            // if (it == nrelfex) p = p0;

            grad_p = esti_grad(p);
            v = v - (step/2) * grad_p;
        }
        h2 = Hamiltonian<NT>(p, v);
        // std::cout << "h2 = " << h2 << std::endl;
        u = ((double)rand()/(double)RAND_MAX);
        if (u < std::exp(h1-h2)) p0 = p;
    }
    return p;
}

template
<
    typename NT,
    typename PointType,
    typename PointList,
    typename RNGType
>
void gaussian_correlation_sampling(CorreSpectra<PointType> &P,
                                    PointList &randPoints,
                                    RNGType &rng,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    const PointType &starting_point,
                                    const NT &step,
                                    unsigned int const& nburns = 0){
    PointType p = starting_point;
    int i = 0;
    while(i < nburns) {
        p = ReHMC(P, p, walkL, rng, step);
        ++i;
    }
    for(i = 0; i < num_points; ++i){
        p = ReHMC(P, p, walkL, rng, step);
        randPoints.push_back(p);
    }   
}

// template
// <
//     typename WalkType
// >
// struct ExponentialCorrelationMatrices
// {
//     template
//     <
//         typename Point,
//         typename NT,
//         typename PointList,
//         typename WalkPolicy,
//         typename RandomNumberGenerator
//     >
//     static void apply(CorreSpectra<Point> const& P,
//                       Point &p,   // a point to start
//                       Point const& c,   // bias function
//                       NT const& T, // temperature/variance
//                       unsigned int const& rnum,
//                       unsigned int const& walk_length,
//                       PointList &randPoints,
//                       WalkPolicy &policy,
//                       RandomNumberGenerator &rng)
//     {
//         Walk walk(P, p, c, T, rng);
//         bool success;
//         for (unsigned int i=0; i<rnum; ++i)
//         {
//             success = walk.template apply(P, p, walk_length, rng);
//             if (!success) {
//                 //return;
//                 throw std::range_error("A generated point is outside polytope");
//             }
//             policy.apply(randPoints, p);
//         }
//     }
// }