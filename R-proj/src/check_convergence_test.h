
#ifndef CHECK_CONVERGENCE_TEST_HPP
#define CHECK_CONVERGENCE_TEST_HPP

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>


template <typename NT>
std::pair<NT, NT> get_mean_variance(std::vector<NT>& vec)
{
    NT mean = 0;
    NT M2 = 0;
    NT variance = 0;
    NT delta;

    unsigned int i=0;
    for (auto vecit = vec.begin(); vecit!=vec.end(); vecit++, i++)
    {
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }
    return std::pair<NT, NT> (mean, variance);
}


template <typename VT, typename NT, typename ConvexBody, typename MT, typename RNGType>
std::pair< std::pair<bool,bool>, std::pair<NT, VT> > check_convergence_test(ConvexBody const& P,
                       MT const& randPoints,
                       bool& too_few,
                       NT& ratio,
                       int const& nu,
                       bool const& lastball,
                       NT const& alpha,
                       NT const& lb,
                       NT const& ub,
                       RNGType &rng)
{
    //NT alpha = parameters.alpha;
    std::pair< std::pair<bool,bool>, std::pair<NT, VT> > res;
    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int NN = randPoints.cols();
    int m = NN/nu;
    NT T;
    NT rs;
    NT alpha_check = 0.01;
    size_t countsIn = 0, countsIn_total = 0;
    bool precheck = false;
    

    for (int i=0; i<NN; i++)
    {
        if (P.is_in(randPoints.col(i))==-1){
            countsIn++;
            countsIn_total++;
            if (rng.sample_urdist() < (NT(1) / countsIn_total))
            {
                res.second.second = randPoints.col(i);
            }
        }
        if ((i+1) % m == 0)
        {
            ratios.push_back(NT(countsIn)/m);
            countsIn = 0;
            //counter++;
            /*if (ratios.size() > 2 && !lastball)
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
                    res.first.second = too_few;
                    res.first.first = false;
                    return res;
                } else if (ratio - T > ub){
                    res.first.first = false;
                    return res;
                }
            }*/
        }
    }

    //if (precheck) alpha *= 0.5;
    mv = get_mean_variance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs * (boost::math::quantile(boost::math::complement(dist, alpha))
           / std::sqrt(NT(nu)));

    res.second.first = ratio;

    if (ratio > lb + T)
    {
        if (lastball)
        {
            res.first.first = true;
            res.first.second = too_few;
            return res;
        }
        if (ratio < ub + T)
        {
            res.first.first = true;
            res.first.second = too_few;
            return res;
        }
        res.first.first = false;
        res.first.second = too_few;
        return res;
    }
    too_few = true;
    res.first.first = false;
    res.first.second = too_few;
    return res;
}


#endif
