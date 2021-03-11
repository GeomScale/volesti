// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file


#ifndef QUADRATIC_POLYNOMIAL_SOLVERS_H
#define QUADRATIC_POLYNOMIAL_SOLVERS_H

// The function compute the roots of a quadratic polynomial equation
template <typename NT>
void solve_quadratic_polynomial(NT a, NT b, NT c, NT &x1, NT &x2, bool &real) {

    real = true;

    if (a == NT(0)) {
        x1 = -c / b;
        x2 = x1;
        return;
    }

    NT Delta = b * b - 4.0 * a * c;
    if (Delta < NT(0)) {
        real = false;
        return;
    }

    if (b >= NT(0)){
        x1 = (- b - std::sqrt(Delta)) / (2.0 * a);
        x2 = (2.0 * c) / (- b - std::sqrt(Delta));
    } else {
        x1 = (2.0 * c) / (- b + std::sqrt(Delta));
        x2 = (- b + std::sqrt(Delta)) / (2.0 * a);
    }
}


#endif

