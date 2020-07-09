// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "EnvelopeProblemSOS.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/fmt/ostr.h"
#include "ChebTools/ChebTools.h"

EnvelopeProblemSOS::EnvelopeProblemSOS(unsigned num_variables, unsigned max_degree, HyperRectangle &hyperRectangle_) :
        _n(num_variables), _d(max_degree),
        _hyperRectangle(hyperRectangle_) {
    assert(num_variables == hyperRectangle_.size());
    //FIXME: for now only univariate polynomials.
    assert(_n == 1);

    _logger = spdlog::get("EnvelopeProblemSOS");
    if (_logger == nullptr) {
        std::vector<spdlog::sink_ptr> sinks;
        sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
        sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/logfile.txt"));
        _logger = std::make_shared<spdlog::logger>("EnvelopeProblemSOS", begin(sinks), end(sinks));
        _logger->set_level(spdlog::level::info);
    }

    _L = _d + 1;
    _U = 2 * _d + 1;

    InterpolantDualSOSBarrier aux_interpolant_barrier(_d);
    std::vector<InterpolantDouble> &chebyshev_points = aux_interpolant_barrier.get_basis();
    InterpolantDouble *cheb_ptr = &chebyshev_points[0];
    InterpolantVector cheb_vec = Eigen::Map<InterpolantVector>(cheb_ptr, chebyshev_points.size());

//    cxxtimer::Timer interp_basis_timer;
//    interp_basis_timer.start();
//    _basis_polynomials.resize(_U);
//    for (unsigned l = 0; l < _basis_polynomials.size(); ++l) {
//        _basis_polynomials[l] = InterpolantVector::Zero(_U);
//    }
//    InterpolantVector aux_vec(_U);
//    for (unsigned i = 0; i < _U; ++i) {
//        _logger->debug("Construct {}-th basis element.", i);
//        InterpolantVector &poly_i = _basis_polynomials[i];
//        poly_i(0) = 1;
//        InterpolantDouble denom = 1.;
//        for (unsigned j = 0; j < _U; ++j) {
//            if (i != j) {
//                denom *= chebyshev_points[i] - chebyshev_points[j];
//                aux_vec = -chebyshev_points[j] * poly_i;
////                aux_vec(1) = 1.;
//                InterpolantVector shift(_U + 1);
//                shift << 0, poly_i;
//                poly_i = shift.segment(0, _U) + aux_vec;
////                poly_i = prod_sos(poly_i, fac_j);
//            }
//        }
//        poly_i /= denom;
////        _logger->info("{}-th self constructed interpolant polynomial is: \n{}", i, poly_i.transpose());
//
////        _basis_polynomials[i] = poly_i;
//    }
//    interp_basis_timer.stop();
//    _logger->info("Finished construction in {} seconds.",
//                  interp_basis_timer.count<std::chrono::milliseconds>() / 1000.);


    for (unsigned k = 0; k < _basis_polynomials.size(); ++k) {
        _logger->debug("The {}-th polynomial is:", k);
        for (unsigned i = 0; i < _basis_polynomials[k].size(); ++i) {
            _logger->debug("{}", _basis_polynomials[k]);
        }
    }

    _logger->info("Construct objectives vector...");

    //Old way of computing the objective
//    _objectives_vector.resize(_U);
//    for (unsigned i = 0; i < _U; ++i) {
//        InterpolantVector &poly = _basis_polynomials[i];
//        InterpolantDouble obj = 0;
//        for (unsigned j = 0; j < _U; ++j) {
//            InterpolantDouble upper_bound_term = poly(j) * pow(_hyperRectangle[0].second, j + 1) / (j + 1);
//            InterpolantDouble lower_bound_term = poly(j) * pow(_hyperRectangle[0].first, j + 1) / (j + 1);
//            obj += upper_bound_term - lower_bound_term;
//        }
//        _objectives_vector(i) = -obj;
//    }

    //Faster, cleverer Clenshaw-Curtis algorithm

    get_clenshaw_curtis_integrals();

}

void EnvelopeProblemSOS::get_clenshaw_curtis_integrals() {

    //TODO: Check if _U even is necessary and how to fix for N odd

    //TODO: Speedup by only computing lower / upper half

    Matrix D(_L, _L);
    for (int k = 0; k < _L; ++k) {
        for (int n = 0; n < _L; ++n) {
            const double scale_fac = (n == 0 or n == _L - 1) ? .5 : 1;
            double cos_kn = cos(k * n * EIGEN_PI / (_L - 1));
            D(k, n) = cos_kn * scale_fac;
        }
    }
    D /= (_L - 1);

    Vector FourierCoeff(_L);
    FourierCoeff(0) = 1;
    FourierCoeff(_L - 1) = 1. / (1 - int((_U - 1) * (_U - 1)));

    for (int m = 1; m < _L - 1; ++m) {
        FourierCoeff(m) = 2. / (1 - 2 * m * 2 * m);
    }

    Vector ClenshawCurtisWeights(_U);
    ClenshawCurtisWeights.segment(0, _L) = D.transpose() * FourierCoeff;
    ClenshawCurtisWeights.segment(_L, _L - 1) = ClenshawCurtisWeights.segment(0, _L - 1).reverse();

    ClenshawCurtisWeights(_L - 1) *= 2;

    IPMDouble dummy_double;
//    InterpolantDouble obj_norm = _objectives_vector.norm();
//    ClenshawCurtisWeights *= InterpolantDoubletoIPMDouble(obj_norm, dummy_double) / ClenshawCurtisWeights.norm();

//    Vector obj_double = InterpolantVectortoVector(_objectives_vector, ClenshawCurtisWeights);

//    _logger->info("Objectives difference is {}", (ClenshawCurtisWeights + obj_double).norm());

    _objectives_vector.resize(_U);
    for (int i = 0; i < _objectives_vector.rows(); i++) {
        _objectives_vector(i) = -ClenshawCurtisWeights(i);
    }
}

void EnvelopeProblemSOS::add_polynomial(InterpolantVector &polynomial) {

    InterpolantMatrix Q = get_transformation_matrix();
    _logger->info("Transformation matrix has norm {}", Q.norm());

    if (_input_in_interpolant_basis) {
        //Generate Random polynomial(Might help stability, have to check)
//        polynomial = InterpolantVector::Random(polynomial.rows());
        _polynomials_bounds.push_back(polynomial);
        return;
    }

    _logger->info("Invert transformation matrix to add polynomial ...");
    cxxtimer::Timer invert_timer;
    invert_timer.start();
    InterpolantMatrix Q_inv = Q.inverse();
    invert_timer.stop();
    _logger->info("Inversion took {} seconds.", invert_timer.count<std::chrono::milliseconds>() / 1000.);
    auto inv_error = (Q * Q_inv - InterpolantMatrix::Identity(Q.rows(), Q.cols())).norm();
    _logger->info("Inversion error is {}", inv_error);
    //FIXME: sys solve below does not work. Figure out why.
//    auto inv_sol = Q_inv * polynomial;
    cxxtimer::Timer sys_solve_timer;
    sys_solve_timer.start();
    InterpolantVector inv_sol = Q.colPivHouseholderQr().solve(polynomial);
    sys_solve_timer.stop();
    _logger->info("Solving system took {} seconds.",
                  sys_solve_timer.count<std::chrono::milliseconds>() / 1000.);


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

    for (unsigned poly_idx = 0; poly_idx < NUM_POLYNOMIALS - 1; ++poly_idx) {
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

    for (unsigned poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
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
    InterpolantVector solution(_U);
    if (not _input_in_interpolant_basis) {
        InterpolantMatrix Q = get_transformation_matrix();
        InterpolantVector solution = Q * sol_vec;
    } else {
        solution = sol_vec;
    }
//    std::cout << "Solution polynomial is: " << solution.transpose() << std::endl;
}

void EnvelopeProblemSOS::plot_polynomials_and_solution(const Solution &sol) {

    std::cout << "Create picture of solution. Saved in plot.png..." << std::endl;
    int num_points = 1000;
    assert(num_points > 1);
    std::vector<double> x(num_points);
    assert(_hyperRectangle.size() == 1);

    IPMDouble x_min = _hyperRectangle[0].first;
    IPMDouble x_max = _hyperRectangle[0].second;

    //Cut off boundary because of oscillatory behaviour
    const IPM_DOUBLE delta_x = x_max - x_min;
    x_min += .3 * delta_x;
    x_max -= .3 * delta_x;

    for (int j = 0; j < num_points; ++j) {
        IPMDouble d = x_min + j * (x_max - x_min) / (num_points - 1);
        Double dummy_D;
        x[j] = InterpolantDoubletoIPMDouble(d, dummy_D);
    }

    std::vector<InterpolantVector> poly_plots;
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
    for (unsigned poly_idx = 0; poly_idx < plots.size(); poly_idx++) {
        plots[poly_idx].resize(num_points);
        InterpolantVector poly_in_orig_basis(_U);
        if (_input_in_interpolant_basis) {
            poly_in_orig_basis = poly_plots[poly_idx];
        } else {
            poly_in_orig_basis = Q_interp * poly_plots[poly_idx];
        }
        for (int i = 0; i < num_points; ++i) {
            //Evaluation of vector polynomial at a certain point. Extract method.
            assert(poly_in_orig_basis.size() > 0);
            InterpolantDouble eval = poly_in_orig_basis[0];
            for (int j = 1; j < poly_in_orig_basis.size(); ++j) {
                eval += poly_in_orig_basis[j] * pow(x[i], j);
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

    for (unsigned plt_idx = 0; plt_idx < poly_plots.size() - 1; ++plt_idx) {
        for (unsigned i = 0; i < plots[plt_idx].size(); ++i) {
            y_min = std::min(y_min, plots[plt_idx][i]);
            y_max = std::max(y_max, plots[plt_idx][i]);
        }
    }

    for (unsigned k = 0; k < envelope_plot.size(); ++k) {
        offset_envelope[k] = envelope_plot[k] - (y_max - y_min) / 100.;
    }

    auto y_bound_offset = (y_max - y_min) / 50;
    plt::ylim(y_min - y_bound_offset, y_max + y_bound_offset);
    for (unsigned p_idx = 0; p_idx < poly_plots.size() - 1; ++p_idx) {
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
