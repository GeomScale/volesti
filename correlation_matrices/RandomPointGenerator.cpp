template
<
        typename WalkTypePolicy,
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename NT,
        typename Point
>
void gaussian_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns)
{

    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    PushBackWalkPolicy push_back_policy;

    Point p = starting_point;

    typedef GaussianRandomPointGenerator <walk> RandomPointGenerator;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, p, a, nburns, walk_len, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, p, a, rnum, walk_len, randPoints,
                                push_back_policy, rng);


}

template
<
    typename Walk
>
struct GaussianRandomPointGenerator
{
    template
    <
            typename Polytope,
            typename Point,
            typename NT,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, a_i, rng);

        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};