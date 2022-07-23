template
<
    typename WalkTypePolicy,
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