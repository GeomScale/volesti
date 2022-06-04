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