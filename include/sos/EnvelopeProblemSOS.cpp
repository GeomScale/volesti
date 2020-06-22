#include "EnvelopeProblemSOS.h"

EnvelopeProblemSOS::EnvelopeProblemSOS(int num_variables, int max_degree, HyperRectangle &hyperRectangle_) :
        _n(num_variables), _d(max_degree),
        _hyperRectangle(hyperRectangle_) {
    assert(num_variables == hyperRectangle_.size());
    //FIXME: for now only univariate polynomials.
    assert(_n == 1);

    _L = _d + 1;
    _U = 2 * _d + 1;

    InterpolantDualSOSBarrier aux_interpolant_barrier(_d);
    auto chebyshev_points = aux_interpolant_barrier.get_basis();

    std::cout << "The chebyshev points are: ";

    for (auto cheb : chebyshev_points) {
        std::cout << " " << cheb;
    }

    std::cout << std::endl;

    for (int i = 0; i < _U; ++i) {
        PolynomialSOS p = Vector::Zero(_U);
        p(0) = 1;
        Double denom = 1.0;
        for (int j = 0; j < _U; ++j) {
            if (i != j) {
                denom *= chebyshev_points[i] - chebyshev_points[j];
                PolynomialSOS fac_j = Vector::Zero(_U);
                fac_j(0) = -chebyshev_points[j];
                fac_j(1) = 1;
                p = prod_sos(p, fac_j);
            }
        }
        p /= denom;
        _basis_polynomials.push_back(p);
    }

    for (int k = 0; k < _basis_polynomials.size(); ++k) {
        std::cout << "The " << k << "-th polynomial is:";
        for (int i = 0; i < _basis_polynomials[k].size(); ++i) {
            std::cout << " " << _basis_polynomials[k][i];
        }
        std::cout << std::endl;
    }

    _objectives_vector.resize(_U);
    for (int i = 0; i < _U; ++i) {
        PolynomialSOS &poly = _basis_polynomials[i];
        Double obj = 0;
        for (int j = 0; j < _U; ++j) {
            Double upper_bound_term = poly(j) * pow(_hyperRectangle[0].second, j + 1) / (j + 1);
            Double lower_bound_term = poly(j) * pow(_hyperRectangle[0].first, j + 1) / (j + 1);
            obj += upper_bound_term - lower_bound_term;
        }
        _objectives_vector[i] = -obj;
    }
}

void EnvelopeProblemSOS::add_polynomial(PolynomialSOS &polynomial) {
    auto Q = get_transformation_matrix();
    auto Q_inv = Q.inverse();
    std::cout << "Add polynomial with interpolant transformation matrix " << std::endl;
    std::cout << Q << std::endl;
    std::cout << "And inverse matrix " << std::endl;
    std::cout << Q_inv << std::endl;
    _polynomials_bounds.push_back(Q_inv * polynomial);
}

PolynomialSOS EnvelopeProblemSOS::generate_zero_polynomial(){
    return PolynomialSOS::Zero(_U);
}

Instance EnvelopeProblemSOS::construct_SOS_instance() {
    unsigned const NUM_POLYNOMIALS = _polynomials_bounds.size();
    unsigned const VECTOR_LENGTH = _U;

    //Below commented out is the simple and for IPM troublesome formulation

//        Constraints constraints;
//        constraints.c = Vector::Zero((NUM_POLYNOMIALS + 1) * VECTOR_LENGTH);
//        constraints.c.block(0, 0, VECTOR_LENGTH, 1) = _objectives_vector;
//
//        constraints.A = Matrix::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH, (NUM_POLYNOMIALS + 1) * VECTOR_LENGTH);
//        constraints.b = Vector::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH);
//
//        for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
//            PolynomialSOS polynomial = _polynomials_bounds[poly_idx];
//            Matrix poly_block = Vector::Ones(VECTOR_LENGTH).asDiagonal();
//            //corresponds to X variables
//            constraints.A.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, VECTOR_LENGTH) = poly_block;
//            //corresponds to Y_i variables
//            constraints.A.block(poly_idx * VECTOR_LENGTH, (poly_idx + 1) * VECTOR_LENGTH, VECTOR_LENGTH,
//                                VECTOR_LENGTH) = poly_block;
//            constraints.b.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, 1) = polynomial;
//        }
//
//        std::cout << "Original system was: " << std::endl;
//        constraints.print();
//
//
//        //Construct Barrier function
//
//        ProductBarrier *productBarrier = new ProductBarrier;
//        auto X_barrier = new FullSpaceBarrier(VECTOR_LENGTH);
//        productBarrier->add_barrier(X_barrier);
//
//        for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
//            auto sos_barrier = new InterpolantDualSOSBarrier(_d);
//            productBarrier->add_barrier(sos_barrier);
//        }

    Constraints constraints;
    constraints.c = Vector::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH);
    constraints.c.block(0, 0, VECTOR_LENGTH, 1) = -_objectives_vector;

    constraints.A = Matrix::Zero((NUM_POLYNOMIALS - 1) * VECTOR_LENGTH, NUM_POLYNOMIALS * VECTOR_LENGTH);
    constraints.b = Vector::Zero((NUM_POLYNOMIALS - 1) * VECTOR_LENGTH);

    for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS - 1; ++poly_idx) {
        PolynomialSOS polynomial = _polynomials_bounds[poly_idx + 1];
        Matrix poly_block = Vector::Ones(VECTOR_LENGTH).asDiagonal();
        //corresponds to X variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, VECTOR_LENGTH) = -poly_block;
        //corresponds to Y_i variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, (poly_idx + 1) * VECTOR_LENGTH, VECTOR_LENGTH,
                            VECTOR_LENGTH) = poly_block;
        constraints.b.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, 1) = polynomial - _polynomials_bounds[0];
    }

    std::cout << "Original system was: " << std::endl;
    constraints.print();


    //Construct Barrier function

    ProductBarrier *productBarrier = new ProductBarrier;

    for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
        auto sos_barrier = new InterpolantDualSOSBarrier(_d);
        productBarrier->add_barrier(sos_barrier);
    }

    auto dual_constraints = constraints.dual_system();


    Instance instance;
    instance.constraints = dual_constraints;
    instance.barrier = productBarrier;
    instance.constraints.print();

    return instance;
};


void EnvelopeProblemSOS::print_solution(Solution sol) {
    //Note: the solution we are looking for is the RHS for the first polynomial - the SOS
    // polynomial in the first constraint
    assert(not _polynomials_bounds.empty());
    Vector v = _polynomials_bounds[0] - sol.s.segment(0, _objectives_vector.rows());
    auto Q = get_transformation_matrix();
    PolynomialSOS solution = Q * v;
    std::cout << "Solution polynomial is: " << solution.transpose() << std::endl;
}

void EnvelopeProblemSOS::plot_polynomials_and_solution(const Solution &sol) {

    std::cout << "Print solution to plot.png" << std::endl;
    int num_points = 1000;
    assert(num_points > 1);
    std::vector<Double> x(num_points);
    assert(_hyperRectangle.size() == 1);
    Double x_min = _hyperRectangle[0].first;
    Double x_max = _hyperRectangle[0].second;
    for (int j = 0; j < num_points; ++j) {
        x[j] = x_min + j * (x_max - x_min) / (num_points - 1);
    }

    std::vector<PolynomialSOS> poly_plots;
    for (PolynomialSOS &poly : _polynomials_bounds) {
        poly_plots.push_back(poly);
    }

    assert(not _polynomials_bounds.empty());
    Vector v = _polynomials_bounds[0] - sol.s.segment(0, _objectives_vector.rows());

    poly_plots.push_back(v);

    auto Q = get_transformation_matrix();

    std::vector<std::vector<Double> > plots(poly_plots.size());
    for (int poly_idx = 0; poly_idx < plots.size(); poly_idx++) {
        plots[poly_idx].resize(num_points);
        auto poly_in_standard_basis = Q * poly_plots[poly_idx];
        for (int i = 0; i < num_points; ++i) {
            //Evaluation of vector polynomial at a certain point. Extract method.
            assert(poly_in_standard_basis.size() > 0);
            Double eval = poly_in_standard_basis[0];
            for (int j = 1; j < poly_in_standard_basis.size(); ++j) {
                eval += poly_in_standard_basis[j] * pow(x[i], j);
            }
            plots[poly_idx][i] = eval;
        }
    }

    plt::figure_size(1200, 780);
//        plt::xlim(x_min, x_max);
//        plt::ylim(-700, 500);

    for (int p_idx = 0; p_idx < poly_plots.size() - 1; ++p_idx) {
        plt::plot(x, plots[p_idx]);
    }

    std::vector<Double> &envelope_plot = plots[plots.size() - 1];
    std::vector<Double> offset_envelope(envelope_plot.size());

    //Get plotted min and max;

    Double y_min = std::numeric_limits<Double>::max();
    Double y_max = std::numeric_limits<Double>::min();

    for (int plt_idx = 0; plt_idx < poly_plots.size(); ++plt_idx) {
        for (int i = 0; i < plots[plt_idx].size(); ++i) {
            y_min = std::min(y_min, plots[plt_idx][i]);
            y_max = std::max(y_max, plots[plt_idx][i]);
        }
    }

    for (int k = 0; k < envelope_plot.size(); ++k) {
        offset_envelope[k] = envelope_plot[k] - (y_max - y_min) / 100.;
    }

    plt::named_plot("lower envelope", x, offset_envelope);

    // Plot a line whose name will show up as "log(x)" in the legend.
    // Add graph title
    plt::title("Lower envelope");
    // Enable legend.
    plt::legend();
    plt::save("plot");
}

Matrix EnvelopeProblemSOS::get_transformation_matrix() {
    Matrix Q(_basis_polynomials.size(), _basis_polynomials.size());
    for (int i = 0; i < Q.rows(); ++i) {
        for (int j = 0; j < Q.rows(); ++j)
            Q(i, j) = _basis_polynomials[j][i];
    }
    return Q;
}
