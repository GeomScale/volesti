#include "EnvelopeProblemSOS.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/fmt/ostr.h"

EnvelopeProblemSOS::EnvelopeProblemSOS(int num_variables, int max_degree, HyperRectangle &hyperRectangle_) :
        _n(num_variables), _d(max_degree),
        _hyperRectangle(hyperRectangle_) {
    assert(num_variables == hyperRectangle_.size());
    //FIXME: for now only univariate polynomials.
    assert(_n == 1);

    _logger = spdlog::get("EnvelopeProblemSOS");
    if (_logger == nullptr) {
        std::vector<spdlog::sink_ptr> sinks;
        sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
        sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("log.txt"));
        _logger = std::make_shared<spdlog::logger>("EnvelopeProblemSOS", begin(sinks), end(sinks));
        _logger->set_level(spdlog::level::info);
    }

    _L = _d + 1;
    _U = 2 * _d + 1;

    InterpolantDualSOSBarrier aux_interpolant_barrier(_d);
    std::vector<InterpolantDouble> &chebyshev_points = aux_interpolant_barrier.get_basis();

//    _logger->debug("The chebyshev points are: ", chebyshev_points);
    _logger->info("Construct interpolant polynomial basis...");

    cxxtimer::Timer interp_basis_timer;
    interp_basis_timer.start();
    for (int i = 0; i < _U; ++i) {
        InterpolantVector p = InterpolantVector::Zero(_U);
        p(0) = 1;
        InterpolantDouble denom = 1.;
        for (int j = 0; j < _U; ++j) {
            if (i != j) {
                InterpolantVector fac_j = InterpolantVector::Zero(_U);
                denom *= chebyshev_points[i] - chebyshev_points[j];
                fac_j(0) = -chebyshev_points[j];
                fac_j(1) = 1.;
                p = prod_sos(p, fac_j);
            }
        }
        p /= denom;
        _basis_polynomials.push_back(p);
    }
    interp_basis_timer.stop();
    _logger->info("Finished construction in {} seconds.",
            interp_basis_timer.count<std::chrono::milliseconds>()/1000.);

    //Test basis

//    InterpolantMatrix test_matrix(_basis_polynomials.size(), chebyshev_points.size());
//    for (int l = 0; l < _basis_polynomials.size(); ++l) {
//        auto &poly = _basis_polynomials[l];
//        for (int i = 0; i < chebyshev_points.size(); ++i) {
//            auto &p = chebyshev_points[i];
//            InterpolantDouble eval = poly(0);
//            for (int j = 1; j < poly.size(); ++j) {
//                eval += poly(j) * pow(p, j);
//            }
//            test_matrix(l, i) = eval;
//        }
//    }
//
//    std::cout << "Test lagrange basis: " << std::endl << test_matrix << std::endl;

    for (int k = 0; k < _basis_polynomials.size(); ++k) {
        _logger->debug("The {}-th polynomial is:", k);
        for (int i = 0; i < _basis_polynomials[k].size(); ++i) {
            _logger->debug("{}", _basis_polynomials[k]);
//            std::cout << " " << _basis_polynomials[k][i];
        }
//        std::cout << std::endl;
    }

    //Alternative way of computing the bases
//    std::vector<InterpolantVector> alternative_basis;
//    for (int i = 0; i < _U; ++i) {
//        InterpolantVector p = InterpolantVector::Zero(_U);
//        p(0) = 1;
//        std::vector<InterpolantVector> poly_vec;
//        poly_vec.push_back(p);
//        for (int j = 0; j < _U; ++j) {
//            if (i != j) {
//                InterpolantVector fac_j = InterpolantVector::Zero(_U);
//                fac_j(0) = -chebyshev_points[j] / (chebyshev_points[i] - chebyshev_points[j]);
//                fac_j(1) = 1 / (chebyshev_points[i] - chebyshev_points[j]);
//                poly_vec.push_back(fac_j);
//            }
//        }
//        p = prod_sos(poly_vec);
//        alternative_basis.push_back(p);
//    }


//    InterpolantMatrix alternative_test_matrix(alternative_basis.size(), chebyshev_points.size());
//    for (int l = 0; l < alternative_basis.size(); ++l) {
//        auto &poly = alternative_basis[l];
//        for (int i = 0; i < chebyshev_points.size(); ++i) {
//            auto &p = chebyshev_points[i];
//            InterpolantDouble eval = poly(0);
//            for (int j = 1; j < poly.size(); ++j) {
//                eval += poly(j) * pow(p, j);
//            }
//            alternative_test_matrix(l, i) = eval;
//        }
//    }
//
//    std::cout << "Test alternative lagrange basis: " << std::endl << alternative_test_matrix << std::endl;


    _logger->info("Construct objectives vector...");
    _objectives_vector.resize(_U);
    for (int i = 0; i < _U; ++i) {
        InterpolantVector &poly = _basis_polynomials[i];
        InterpolantDouble obj = 0;
        for (int j = 0; j < _U; ++j) {
            InterpolantDouble upper_bound_term = poly(j) * pow(_hyperRectangle[0].second, j + 1) / (j + 1);
            InterpolantDouble lower_bound_term = poly(j) * pow(_hyperRectangle[0].first, j + 1) / (j + 1);
            obj += upper_bound_term - lower_bound_term;
        }
        _objectives_vector[i] = -obj;
    }
    _logger->info("Finished constructing objective vector");
}

void EnvelopeProblemSOS::add_polynomial(InterpolantVector &polynomial) {
    InterpolantMatrix Q = get_transformation_matrix();
    InterpolantMatrix Q_inv = Q.inverse();
    auto inv_error = (Q * Q_inv - InterpolantMatrix::Identity(Q.rows(), Q.cols())).norm();
    std::cout << "Inversion error is " << inv_error << std::endl;
//    std::cout << "Add polynomial with interpolant transformation matrix " << std::endl;
//    std::cout << Q << std::endl;
//    std::cout << "And inverse matrix " << std::endl;
//    std::cout << Q_inv << std::endl;
    //FIXME: sys solve below does not work. Figure out why.
//    auto sol = Q.fullPivLu().solve(polynomial);
//    std::cout << "lin sys solution is: " << sol.transpose() << std::endl;
//    double relative_error = (Q*sol - polynomial).norm() / polynomial.norm(); // norm() is L2 norm
//    std::cout << "The relative error is:\n" << relative_error << std::endl;
    auto inv_sol = Q_inv * polynomial;
//    std::cout << "Solutions differ in norm by " << (sol - inv_sol).norm() << std::endl;
//    std::cout << "The errors are: " << std::endl;

//    std::cout << "Originally: " << (polynomial - Q * inv_sol).norm() << std::endl;
//    std::cout << " And new: " << (polynomial - Q * sol).norm() << std::endl;

    _polynomials_bounds.push_back(inv_sol);
}

InterpolantVector EnvelopeProblemSOS::generate_zero_polynomial() {
    return InterpolantVector::Zero(_U);
}

Instance EnvelopeProblemSOS::construct_SOS_instance() {
    unsigned const NUM_POLYNOMIALS = _polynomials_bounds.size();
    unsigned const VECTOR_LENGTH = _U;

    Constraints constraints;
    constraints.c = Vector::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH);
    constraints.c.block(0, 0, VECTOR_LENGTH, 1) = -InterpolantVectortoVector(_objectives_vector,
                                                                             constraints.c);

    constraints.A = Matrix::Zero((NUM_POLYNOMIALS - 1) * VECTOR_LENGTH, NUM_POLYNOMIALS * VECTOR_LENGTH);
    constraints.b = Vector::Zero((NUM_POLYNOMIALS - 1) * VECTOR_LENGTH);

    for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS - 1; ++poly_idx) {
        PolynomialSOS dummy;
        //dummy Polynomial to infer type.
        PolynomialSOS polynomial = InterpolantVectortoVector(_polynomials_bounds[poly_idx + 1], dummy);
        Matrix poly_block = Vector::Ones(VECTOR_LENGTH).asDiagonal();
        //corresponds to X variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, VECTOR_LENGTH) = -poly_block;
        //corresponds to Y_i variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, (poly_idx + 1) * VECTOR_LENGTH, VECTOR_LENGTH,
                            VECTOR_LENGTH) = poly_block;
        constraints.b.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, 1) =
                polynomial - InterpolantVectortoVector(_polynomials_bounds[0], constraints.b);
    }

    std::cout << "Original SOS instance created." << std::endl;
//    constraints.print();

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

    std::cout << "Dual formulation created." << std::endl;
//    instance.constraints.print();

    return instance;
};


void EnvelopeProblemSOS::print_solution(Solution sol) {
    //Note: the solution we are looking for is the RHS for the first polynomial - the SOS
    // polynomial in the first constraint
    assert(not _polynomials_bounds.empty());
    Vector dummy_vec;
    Matrix dummy_matrix;
    Vector sol_seg = sol.s.segment(0, _objectives_vector.rows());
    InterpolantVector sol_seg_interp(sol_seg.rows());
    for (int i = 0; i < sol_seg.rows(); i++) {
        sol_seg_interp(i) = sol_seg(i);
    }
    InterpolantVector sol_vec = _polynomials_bounds[0] - sol_seg_interp;
    InterpolantMatrix Q = get_transformation_matrix();
    InterpolantVector solution = Q * sol_vec;
//    std::cout << "Solution polynomial is: " << solution.transpose() << std::endl;
}

void EnvelopeProblemSOS::plot_polynomials_and_solution(const Solution &sol) {

    std::cout << "Create solution picture at plot.png..." << std::endl;
    int num_points = 1000;
    assert(num_points > 1);
    std::vector<double> x(num_points);
    assert(_hyperRectangle.size() == 1);
    IPMDouble x_min = _hyperRectangle[0].first;
    IPMDouble x_max = _hyperRectangle[0].second;
    for (int j = 0; j < num_points; ++j) {
        IPMDouble d = x_min + j * (x_max - x_min) / (num_points - 1);
        Double dummy_D;
        x[j] = InterpolantDoubletoIPMDouble(d, dummy_D);
    }

    std::vector<InterpolantVector> poly_plots;
    IPMDouble dummy_ipm_d;
    for (InterpolantVector &poly : _polynomials_bounds) {
        Vector poly_double(poly.rows());
        poly_plots.push_back(poly);
    }

    assert(not _polynomials_bounds.empty());

    Vector seg_vec = sol.s.segment(0, _objectives_vector.rows());
    InterpolantVector seg_interp_vec(seg_vec.rows());
    for (int i = 0; i < seg_vec.rows(); i++) {
        seg_interp_vec(i) = seg_vec(i);
    }

    InterpolantVector v = _polynomials_bounds[0] - seg_interp_vec;

    poly_plots.push_back(v);

    InterpolantMatrix Q_interp = get_transformation_matrix();
    Matrix dummy_matrix;
//    Matrix Q = InterpolantMatrixToMatrix(Q_interp,dummy_matrix);

    std::vector<std::vector<double> > plots(poly_plots.size());
    for (int poly_idx = 0; poly_idx < plots.size(); poly_idx++) {
        plots[poly_idx].resize(num_points);
        InterpolantVector poly_in_standard_basis = Q_interp * poly_plots[poly_idx];
        for (int i = 0; i < num_points; ++i) {
            //Evaluation of vector polynomial at a certain point. Extract method.
            assert(poly_in_standard_basis.size() > 0);
            InterpolantDouble eval = poly_in_standard_basis[0];
            for (int j = 1; j < poly_in_standard_basis.size(); ++j) {
                eval += poly_in_standard_basis[j] * pow(x[i], j);
            }
            Double dummy_D;
            plots[poly_idx][i] = InterpolantDoubletoIPMDouble(eval, dummy_D);
        }
    }

    plt::figure_size(1200, 780);
//        plt::xlim(x_min, x_max);
//        plt::ylim(-700, 500);



    std::vector<double> &envelope_plot = plots[plots.size() - 1];
    std::vector<double> offset_envelope(envelope_plot.size());

    //Get plotted min and max;

    double y_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::min();

    for (int plt_idx = 0; plt_idx < poly_plots.size() - 1; ++plt_idx) {
        for (int i = 0; i < plots[plt_idx].size(); ++i) {
            y_min = std::min(y_min, plots[plt_idx][i]);
            y_max = std::max(y_max, plots[plt_idx][i]);
        }
    }

    for (int k = 0; k < envelope_plot.size(); ++k) {
        offset_envelope[k] = envelope_plot[k] - (y_max - y_min) / 100.;
    }

    auto y_bound_offset = (y_max - y_min) / 50;
    plt::ylim(y_min - y_bound_offset, y_max + y_bound_offset);
    for (int p_idx = 0; p_idx < poly_plots.size() - 1; ++p_idx) {
        plt::plot(x, plots[p_idx]);
    }
    plt::named_plot("lower envelope", x, offset_envelope);

    // Plot a line whose name will show up as "log(x)" in the legend.
    // Add graph title
    plt::title("Lower envelope");
    // Enable legend.
    plt::legend();
    plt::save("plot");

    std::cout << "Done." << std::endl;
}

InterpolantMatrix EnvelopeProblemSOS::get_transformation_matrix() {
    InterpolantMatrix Q(_basis_polynomials.size(), _basis_polynomials.size());
    for (int i = 0; i < Q.rows(); ++i) {
        for (int j = 0; j < Q.rows(); ++j)
            Q(i, j) = _basis_polynomials[j][i];
    }
    return Q;
}
