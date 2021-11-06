
#ifndef CHECK_CONVERGENCE_HPP
#define CHECK_CONVERGENCE_HPP

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>


template <typename ConvexBody, typename PointList, typename NT>
bool check_convergence(ConvexBody const& P,
                       MT const& randPoints,
                       bool& too_few,
                       NT& ratio,
                       int const& nu,
                       bool const& precheck,
                       bool const& lastball,
                       NT const& alpha,
                       NT const& lb,
                       NT const& ub)
{
    //NT alpha = parameters.alpha;
    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.cols()/nu;
    NT T;
    NT rs;
    NT alpha_check = 0.01;
    size_t countsIn = 0;

    for (int i=0; i < randPoints.cols(); i++)
    {
        if (P.is_in(randPoints.col(i))==-1) countsIn++;
        if (i % m == 0)
        {
            ratios.push_back(NT(countsIn)/m);
            countsIn = 0;
            if (ratios.size()>1 && precheck)
            {
                boost::math::students_t dist(ratios.size() - 1);
                mv = get_mean_variance(ratios);
                ratio = mv.first;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile
                            (boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < lb)
                {
                    too_few = true;
                    return false;
                } else if (ratio - T > ub) return false;
            }
        }
    }

    if (precheck) alpha *= 0.5;
    mv = get_mean_variance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs * (boost::math::quantile(boost::math::complement(dist, alpha))
           / std::sqrt(NT(nu)));
    if (ratio > lb + T)
    {
        if (lastball) return true;
        if ((precheck && ratio < ub - T)
        || (!precheck && ratio < ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;
}


#endif
