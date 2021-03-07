


#ifndef QUADRATIC_POLYNOMIAL_ROOTS_H
#define QUADRATIC_POLYNOMIAL_ROOTS_H

template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template <typename NT>
void solve_qudratic_polynomial(NT a, NT b, NT c, NT &x1, NT &x2, bool &real)
{
    a = NT(1);
    b = b / a;
    c = c / a;
    real = true;

    NT sqrt_c = std::sqrt(std::abs(c));
    NT alpha = sgn(b) * sqrt_c;
    NT e = NT(sgn(c));
    NT beta = std::abs(b) / (NT(2)*sqrt_c);

    if (e > NT(0) && beta < NT(1))
    {
        real = false;
        return;
    }

    if (e < NT(0))
    {
        x1 = beta + std::sqrt(beta*beta + NT(1));
        x2 = -NT(1)/x1;
    } else 
    {
        x1 = beta + std::sqrt((beta + NT(1)) * (beta - NT(1)));
        x2 = NT(1) / x1;
    }
    x1 = -x1 * alpha;
    x2 = -x2 * alpha;
}



#endif

