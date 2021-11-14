
#ifndef SAMPLING_SEGMENT_HPP
#define SAMPLING_SEGMENT_HPP

#include <cmath>

template <typename NT, typename RandomNumberGenerator>
NT sample_sigma_gaussian_segment(NT const& t1, NT const& t2, NT const& a, NT const& b, NT const& c, 
                                 NT const& d, NT const& e, unsigned int const& W, NT const &k, 
                                 RandomNumberGenerator &rng)
{

    NT t = (t1+t2)*0.5;
    NT L = (0.33*(t2-t1))*0.5;
    NT t_cos, t_sin, t_next, py, l;

    t_cos = cos(t);
    t_sin = sin(t);

    NT px = a*(t_cos*t_cos) + b*(t_sin*t_sin) + c*t_cos*t_sin + d*t_sin + e*t_cos;

    for (int i=0; i<W; i++)
    {
        l = rng.sample_urdist();
        t_next = l * (t + L) + (1 - l) * (t - L);

        if (t_next>t2 || t_next<t1)
        {
            continue;
        }

        t_cos = cos(t_next);
        t_sin = sin(t_next);

        py = a*(t_cos*t_cos) + b*(t_sin*t_sin) + c*t_cos*t_sin + d*t_sin + e*t_cos;

        if (log(rng.sample_urdist()) < -k*(py-px))
        {
            t = t_next;
            px = py;
        }
    }

    return t;
}






template <typename NT, typename RandomNumberGenerator>
std::pair<NT, bool> sample_sigma_gaussian_segment_with_check(NT const& t1, NT const& t2, NT const& a, NT const& b, NT const& c, 
                                            NT const& d, NT const& e, unsigned int const& W, NT const &k, 
                                            RandomNumberGenerator &rng)
{

    NT t = (t1+t2)*0.5;
    NT L = (0.33*(t2-t1))*0.5;
    NT t_cos, t_sin, t_next, py, l;
    bool got_outside = false;

    t_cos = cos(t);
    t_sin = sin(t);

    NT px = a*(t_cos*t_cos) + b*(t_sin*t_sin) + c*t_cos*t_sin + d*t_sin + e*t_cos;

    for (int i=0; i<W; i++)
    {
        l = rng.sample_urdist();
        t_next = l * (t + L) + (1 - l) * (t - L);

        

        if (t_next>t2 || t_next<t1)
        {
            if (!got_outside)
            {
                t_cos = cos(t_next);
                t_sin = sin(t_next);

                py = a*(t_cos*t_cos) + b*(t_sin*t_sin) + c*t_cos*t_sin + d*t_sin + e*t_cos;
            
                if (log(rng.sample_urdist()) < -k*(py-px))
                {
                    got_outside = true;
                }
            }
            continue;
        }

        t_cos = cos(t_next);
        t_sin = sin(t_next);

        py = a*(t_cos*t_cos) + b*(t_sin*t_sin) + c*t_cos*t_sin + d*t_sin + e*t_cos;

        if (log(rng.sample_urdist()) < -k*(py-px))
        {
            t = t_next;
            px = py;
        }
    }

    return std::pair<NT, bool>(t, got_outside);
}


#endif
