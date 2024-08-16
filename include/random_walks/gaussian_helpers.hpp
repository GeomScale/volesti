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

#endif // GAUSSIAN_HELPERS_HPP
