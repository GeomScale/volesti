// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Based on https://papers.nips.cc/paper/8483-the-randomized-midpoint-method-for-log-concave-sampling.pdf

#ifndef ODE_SOLVERS_RANDOMIZED_MIDPOINT_HPP
#define ODE_SOLVERS_RANDOMIZED_MIDPOINT_HPP

template <
typename Point,
typename NT,
typename Polytope,
typename func,
typename RandomNumberGenerator
>
struct RandomizedMipointSDESolver {

    typedef std::vector<Point> pts;

    typedef std::vector<Polytope*> bounds;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned int dim;

    VT Ar, Av;

    NT eta;
    NT t;
    NT u;

    func F;
    bounds Ks;

    // Contains the sub-states
    pts xs;
    pts xs_prev;
    Point y, w, z;

    Point W1, W2, W3;
    VT Y, Z;
    MT C, L, D, M;

    std::pair<NT, int> pbpair;

    RandomizedMipointSDESolver(NT initial_time, NT step, pts initial_state, func oracle, bounds boundaries, NT u_=NT(1)) :
    eta(step), t(initial_time), u(u_), F(oracle), Ks(boundaries), xs(initial_state) {
        dim = xs[0].dimension();
        Y.resize(2 * dim);
        Z.resize(2 * dim);
        C.resize(2 * dim, 2 * dim);
        L.resize(2 * dim, 2 * dim);
        D.resize(2 * dim, 2 * dim);
        M.resize(2 * dim, 2 * dim);
        W1 = Point(dim);
        W2 = Point(dim);
        W3 = Point(dim);

        for (unsigned int i = 0; i < 2 * dim; i++) {
            Y(i) = NT(0);
            Z(i) = NT(0);
            for (unsigned int j = 0; j < 2 * dim; j++) {
                C(i, j) = NT(0);
                D(i, j) = NT(0);
            }
        }
    };

    void step(RandomNumberGenerator &rng) {
        xs_prev = xs;
        unsigned int x_index, v_index, it;
        NT a;
        t += eta;
        for (unsigned int i = 1; i < xs.size(); i += 2) {

            a = rng.sample_urdist();

            x_index = i - 1;
            v_index = i;

            calculate_Ws(a, rng);

            z = xs_prev[x_index];
            z = z + (0.5 * (1 - exp(- 2 * a * eta))) * xs_prev[v_index];
            z = z - (0.5 * u * (a * eta - 0.5 * (1 - exp(-2 * a * eta)))) * F(v_index, xs_prev, t);
            z = z + sqrt(u) * W1;

            w = xs_prev[x_index];
            xs[x_index] = xs_prev[x_index];
            y = (-1.0) * xs_prev[x_index];
            xs_prev[x_index] = z;
            xs[x_index] = xs[x_index] + (0.5 * (1 - exp(-2 * eta))) * xs[v_index];
            xs[x_index] = xs[x_index] + (0.5 * u * eta * (1 - exp(-2 * (eta - a * eta)))) * F(v_index, xs_prev, t);
            xs[x_index] = xs[x_index] + sqrt(u) * W2;

            xs[v_index] = exp(- 2 * eta) * xs_prev[v_index];
            xs[v_index] = xs[v_index] + u * eta * exp(-2 * (eta - a * eta)) * F(v_index, xs_prev, t);
            xs[v_index] = xs[v_index] + (2 * sqrt(u)) * W3;

            xs_prev[x_index] = w;

            y = y + xs[x_index];

            if (Ks[x_index] == NULL) {
                xs[x_index] = xs_prev[x_index] + y;
            }
            else {
                // Find intersection (assuming a line trajectory) between x and y
                do {

                    pbpair = Ks[x_index]->line_positive_intersect(xs_prev[x_index], y, Ar, Av);

                    if (pbpair.first >= 0 && pbpair.first <= 1) {
                        xs_prev[x_index] += (pbpair.first * 0.95) * y;
                        Ks[x_index]->compute_reflection(y, xs_prev[x_index], pbpair.second);
                        xs[x_index] = xs_prev[x_index] + y;

                        // Reflect velocity
                        Ks[x_index]->compute_reflection(xs[v_index], xs[x_index], pbpair.second);
                    }
                    else {
                        xs[x_index] = xs_prev[x_index] + y;
                    }
                } while (!Ks[x_index]->is_in(xs[x_index]));
            }

        }
    }

    void print_state() {
        for (int j = 0; j < xs.size(); j ++) {
            for (unsigned int i = 0; i < xs[j].dimension(); i++) {
                std::cout << xs[j][i] << " ";
            }
        }
        std::cout << std::endl;
    }

    void steps(int num_steps, RandomNumberGenerator &rng) {
        for (int i = 0; i < num_steps; i++) step(rng);
    }

    Point get_state(int index) {
        return xs[index];
    }

    void set_state(int index, Point p) {
        xs[index] = p;
    }

    void calculate_Ws(NT &a, RandomNumberGenerator &rng) {
        // Initialize matrices to zero
        Y = 0 * Y;
        Z = 0 * Z;
        C = 0 * C;
        D = 0 * D;

        // Helper variables
        NT temp_y, temp_z, h1, h2, g1, g2;

        // Variance of variable G1, G2
        temp_y = 0.25 * (exp(4 * a * eta) - 1);
        temp_z = 0.25 * (exp(4 * eta) - exp(4 * a * eta));
        for (unsigned int i = 0; i < dim; i++) {
            C(i, i) = temp_y;
            D(i, i) = temp_z;
        }

        // Variance of H1, H2
        temp_y = a * eta;
        temp_z = (eta - a * eta);
        for (unsigned int i = dim; i < 2 * dim; i++) {
            C(i, i) = temp_y;
            D(i, i) = temp_z;
        }

        // Covariances of Hi, Gi
        temp_y = 0.5 * (exp(2 * a * eta) - 1);
        temp_z = 0.5 * (exp(2 * eta) - exp(2 * a * eta));
        for (unsigned int i = 0; i < dim; i++) {
            C(i, i + dim) = temp_y;
            C(i + dim, i) = temp_y;
            D(i, i + dim) = temp_z;
            D(i + dim, i) = temp_z;
        }

        // Cholesky Decomposition
        Eigen::LLT<MT> lltofC(C);
        L = lltofC.matrixL();
        Eigen::LLT<MT> lltofD(D);
        M = lltofD.matrixL();

        // Normal Vectors
        for (unsigned int i = 0; i < 2 * dim; i++) {
            Y(i) = rng.sample_ndist();
            Z(i) = rng.sample_ndist();
        }

        // Transformed vectors
        Y = L * Y;
        Z = M * Z;

        // Calculate Brownian Integrals W1, W2, W3 (Appendix A)
        for (int i = 0; i < dim; i++) {
            h1 = Y(i + dim);
            g1 = Y(i);
            h2 = Z(i + dim);
            g2 = Z(i);
            W1.set_coord(i, h1 - exp(- 2 * a * eta) * g1);
            W2.set_coord(i, h1 + h2 - exp(- 2 * a * eta) * (g1 + g2));
            W3.set_coord(i, exp(- 2 * eta) * (g1 + g2));
        }

    }

};


#endif
