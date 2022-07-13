template<typename PointType,
    typename PointList,
    typename RandomNumberGenerator, 
    typename WalkTypePolicy>
MT uniform_correlation_sampling(PointList &randPoints,
                   CorreSpectra P,
                   RandomNumberGenerator &rng,
                   WalkTypePolicy &WalkType,
                   const unsigned int &walk_len,
                   const unsigned int &num_points,
                   unsigned int const& nburns){
    typedef typename WalkTypePolicy::template Walk <CorreSpectra,
                                                    RandomNumberGenerator> walk;
    typedef RandomPointGenerator <walk> RandomPointGenerator;

    PushBackWalkPolicy push_back_policy;
    RandomNumberGenerator rng(P.dimension());
    Point p(P.dimension());
    P.set_interior_point(p);
    
    Generator::apply(HP, p, num_points, walk_len,
                     randPoints, push_back_policy, rng);
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
                                    push_back_policy, rng, WalkType.param);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, num_points, walk_len, randPoints,
                                push_back_policy, rng, WalkType.param);
}