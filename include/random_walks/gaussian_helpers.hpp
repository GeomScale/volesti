#ifndef GAUSSIAN_HELPERS_HPP
#define GAUSSIAN_HELPERS_HPP

#define EXP_CHORD_TOLERENCE 0.00000001

// evaluate the pdf of point p
template <typename Point, typename NT>
NT eval_exp(Point const& p, NT const& a)
{
    return std::exp(-a * p.squared_length());
}
template <typename Point, typename NT, typename MT>
NT eval_exp(Point const& p, MT const& inv_covariance_matrix, NT const& a_next, NT const& a_curr)
{
    Eigen::Matrix<NT, Eigen::Dynamic, 1> dist_vector = p.getCoefficients();
    NT mahalanobis_dist = dist_vector.transpose() * inv_covariance_matrix * dist_vector;
    NT log_ratio = (a_curr - a_next) * mahalanobis_dist;
    return std::exp(log_ratio);
}

template <typename Point, typename NT>
NT get_max(Point const& l, Point const& u, NT const& a_i)
{
    NT res;
    Point a = -1.0 * l;
    Point bef = u - l;
    Point b = (1.0 / std::sqrt((bef).squared_length())) * bef;
    Point z = (a.dot(b) * b) + l;
    NT low_bd = (l[0] - z[0]) / b[0], up_bd = (u[0] - z[0]) / b[0];
    if (low_bd * up_bd > 0)
    {
        //if(std::signbit(low_bd)==std::signbit(up_bd)){
        res = std::max(eval_exp(u, a_i), eval_exp(l, a_i));
    }
    else
    {
        res = eval_exp(z, a_i);
    }

    return res;
}


template <typename NT>
NT get_max_coord(NT const& l, NT const& u, NT const& a_i)
{
    const NT zero = 0;
    if (l < zero && u > zero)
    {
        return NT(1);
    }
    return std::max(std::exp(-a_i * l * l), std::exp(-a_i * u * u));
}


// Pick a point from the distribution exp(-a_i x^2) truncated to [l, u], where
// the interval lies on one side of the mode, i.e. l >= 0 or u <= 0.
// Rejection sampling from a translated exponential proposal for long intervals,
// and from a bounding rectangle for short ones, keeps the acceptance probability
// above a constant however far the interval is from the mode, see
// C. P. Robert, "Simulation of truncated normal variables", Statistics and
// Computing 5, 1995.
template <typename NT, typename RandomNumberGenerator>
NT chord_random_point_generator_exp_tail(NT const& l,
                                         NT const& u,
                                         NT const& a_i,
                                         RandomNumberGenerator& rng)
{
    // reduce to the standard gaussian truncated to [lo, up] with 0 <= lo
    const NT sd = 1.0 / std::sqrt(2.0 * a_i);
    const NT sign = (l >= 0) ? NT(1) : NT(-1);
    const NT lo = ((l >= 0) ? l : -u) / sd;
    const NT up = ((l >= 0) ? u : -l) / sd;

    // the optimal rate of the exponential proposal
    const NT alpha = (lo + std::hypot(lo, NT(2))) / 2.0;
    NT x;
    if ((up - lo) * alpha >= 1.0)
    {
        while (true)
        {
            x = lo - std::log(1.0 - rng.sample_urdist()) / alpha;
            if (x <= up &&
                rng.sample_urdist() <= std::exp(-0.5 * (x - alpha) * (x - alpha)))
            {
                break;
            }
        }
    }
    else
    {
        // the density is maximized at lo
        while (true)
        {
            x = lo + (up - lo) * rng.sample_urdist();
            if (rng.sample_urdist() <= std::exp(-0.5 * (x - lo) * (x + lo)))
            {
                break;
            }
        }
    }
    return sign * sd * x;
}

#endif // GAUSSIAN_HELPERS_HPP
