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
#include "barriers/InterpolantDualSOSBarrier.h"
#include "barriers/ProductBarrier.h"
#include "barriers/SumBarrier.h"
#include <boost/math/special_functions/binomial.hpp>
#include "../../external/Padua/padua.h"

//TODO: use precompiled versions of padua.h and line_fekete_rule.h

template <typename IPMDouble>
EnvelopeProblemSOS<IPMDouble>::EnvelopeProblemSOS(unsigned num_variables, unsigned max_degree, HyperRectangle &hyperRectangle_) :
        _n(num_variables), _d(max_degree),
        _hyperRectangle(hyperRectangle_) {
    assert(num_variables == hyperRectangle_.size());
    //FIXME: for now only univariate polynomials.
    assert(_n == 1);
    initialize_loggers();
    initialize_problem();
}

template <typename IPMDouble>
EnvelopeProblemSOS<IPMDouble>::EnvelopeProblemSOS(std::string  instance_file, std::string config_json) {
    initialize_loggers();

    pt::read_json(instance_file, _instance);
    _n = _instance.get<int>("num_variables");
    _d = _instance.get<int>("max_degree");


    _logger->info("Configure Envelope Problem with file {}", config_json);
    pt::read_json(config_json, _config);

    if(_config.find("use_weighted_polynomials") != _config.not_found()){
        _use_weighted_polynomials = _config.get<bool>("use_weighted_polynomials");
    }

    if(_config.find("input_in_interpolant_basis") != _config.not_found()){
        _input_in_interpolant_basis = _config.get<bool>("input_in_interpolant_basis");
    }

    HyperRectangle hyperRectangle;
    //Note: Keep interval bounds for now.
    IPMDouble const interval_lower_bound = -1.;
    IPMDouble const interval_upper_bound = 1.;
    hyperRectangle.push_back(std::pair<IPMDouble, IPMDouble>(interval_lower_bound,
                                                             interval_upper_bound));

    _hyperRectangle = hyperRectangle;

    initialize_problem();

    for (pt::ptree::value_type &poly : _instance.get_child("polynomials"))
    {
        unsigned idx = 0;
        InterpolantVector sos_poly = generate_zero_polynomial();
        for (pt::ptree::value_type &cell : poly.second)
        {
            sos_poly[idx++] = cell.second.get_value<IPMDouble>();
        }
        add_polynomial(sos_poly);
        _logger->info("Polynomial added");
    }

//    while (std::getline(instance_file, line)) {
//        std::istringstream poly_stream(line);
//        InterpolantVector sos_poly = generate_zero_polynomial();
//        IPMDouble val;
//        unsigned idx = 0;
//        while (poly_stream >> val) {
//            if (idx >= sos_poly.size()) {
//                _logger->error("Input data in wrong format.");
//                exit(1);
//            }
//            sos_poly[idx++] = val;
//        }
//        add_polynomial(sos_poly);
//        _logger->info("Polynomial added.");
//    }
}

template <typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::initialize_problem() {
    _L = static_cast<unsigned>(boost::math::binomial_coefficient<double>(_d + _n, _n));
    _U = static_cast<unsigned>(boost::math::binomial_coefficient<double>(2 * _d + _n, _n));

    _objectives_vector.resize(_U);

    //univariate case
    if(_n == 1) {
        if (not _input_in_interpolant_basis) {
            calculate_basis_polynomials();
        }
        for (unsigned k = 0; k < _basis_polynomials.size(); ++k) {
            _logger->debug("The {}-th polynomial is:", k);
            for (unsigned i = 0; i < _basis_polynomials[k].size(); ++i) {
                _logger->debug("{}", _basis_polynomials[k]);
            }
        }
        _logger->info("Construct objectives vector...");

        //Clenshaw-Curtis algorithm
        compute_clenshaw_curtis_integrals();
    } else if(_n == 2){//Padua points
        double * pd_w = padua::padua_weights(_d + 1);
        for(int i = 0; i < _U; i ++){
            _objectives_vector(i) = -pd_w[i];
        }
    }
    //n >  2
    else {
        //Fekete points
        //TODO: how to do this? I do not understand what Papp and Yildiz are doing here, it seems they
        //use just any objective that guaranteese boundedness.
        //Here we just copy their implementation.

        //TODO: test. Currently not working returning a feasible solution

        InterpolantDualSOSBarrier<IPMDouble> tmp_interp(_d,_n);
        Vector tmp_ones = Vector::Ones(tmp_interp.get_P().cols());
        _objectives_vector = tmp_interp.get_P().colPivHouseholderQr().solve(tmp_ones).template cast<InterpolantDouble>();
        InterpolantVector tmp_zero = InterpolantVector::Zero(_objectives_vector.rows());
        _objectives_vector.cwiseMax(tmp_zero);
//        _objectives_vector = -_objectives_vector;
    }
}


template <typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::initialize_loggers() {
    _logger = spdlog::get("EnvelopeProblemSOS");
    if (_logger == nullptr) {
        std::__1::vector<spdlog::sink_ptr> sinks;
        sinks.push_back(std::__1::make_shared<spdlog::sinks::stdout_color_sink_mt>());
        sinks.push_back(std::__1::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/logfile.txt", true));
        _logger = std::__1::make_shared<spdlog::logger>("EnvelopeProblemSOS", begin(sinks), end(sinks));
        _logger->set_level(spdlog::level::info);
    }
}

template <typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::calculate_basis_polynomials() {

    //Currently only in univariate case.
    assert(_n == 1);

    InterpolantDualSOSBarrier<IPMDouble> aux_interpolant_barrier(_d);
    std::vector<std::vector<InterpolantDouble > > & cheb_vec_tmp = aux_interpolant_barrier.get_basis();
    std::vector<InterpolantDouble> chebyshev_points;
    for(auto vec : cheb_vec_tmp){
        chebyshev_points.push_back(vec[0]);
    }
    InterpolantDouble *cheb_ptr = &chebyshev_points[0];
    InterpolantVector cheb_vec = Eigen::Map<InterpolantVector>(cheb_ptr, chebyshev_points.size());

    _logger->info("Construct transformation matrix");
    cxxtimer::Timer interp_basis_timer;
    interp_basis_timer.start();
    _basis_polynomials.resize(_U);
    for (unsigned l = 0; l < _basis_polynomials.size(); ++l) {
        _basis_polynomials[l] = InterpolantVector::Zero(_U);
    }
    InterpolantVector aux_vec(_U);
    for (unsigned i = 0; i < _U; ++i) {
        _logger->debug("Construct {}-th basis element.", i);
        InterpolantVector &poly_i = _basis_polynomials[i];
        poly_i(0) = 1;
        InterpolantDouble denom = 1.;
        for (unsigned j = 0; j < _U; ++j) {
            if (i != j) {
                denom *= chebyshev_points[i] - chebyshev_points[j];
                aux_vec = -chebyshev_points[j] * poly_i;
                InterpolantVector shift(_U + 1);
                shift << 0, poly_i;
                poly_i = shift.segment(0, _U) + aux_vec;
            }
        }
        poly_i /= denom;
    }
    interp_basis_timer.stop();

    _logger->info("Finished construction in {} seconds.",
                  interp_basis_timer.count<std::chrono::milliseconds>() / 1000.);
}

template <typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::compute_clenshaw_curtis_integrals() {

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

    for (int i = 0; i < _objectives_vector.rows(); i++) {
        _objectives_vector(i) = -ClenshawCurtisWeights(i);
    }
}

template<typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::add_polynomial(InterpolantVector &polynomial) {

    InterpolantMatrix Q = get_transformation_matrix();
    _logger->info("Transformation matrix has norm {}", Q.norm());

    if (_input_in_interpolant_basis) {
        InterpolantDualSOSBarrier<IPMDouble> aux_barrier(_d, _n);
        Matrix P = aux_barrier.get_P();
        InterpolantVector pol = P.template cast<InterpolantDouble>() * polynomial.segment(0, P.cols());
        _polynomials_bounds.push_back(pol);
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
    cxxtimer::Timer sys_solve_timer;
    sys_solve_timer.start();
    InterpolantVector inv_sol = Q.colPivHouseholderQr().solve(polynomial);
    sys_solve_timer.stop();
    _logger->info("Solving system took {} seconds.",
                  sys_solve_timer.count<std::chrono::milliseconds>() / 1000.);
    _polynomials_bounds.push_back(inv_sol);
}

template <typename IPMDouble>
InterpolantVector EnvelopeProblemSOS<IPMDouble>::generate_zero_polynomial() {
    return InterpolantVector::Zero(_U);
}

//TODO: it would be better to introduce barrier functions only in the IPM, not when constructing the instance.
// Although we already have to declare cone membership here.

template<typename IPMDouble>
Instance<IPMDouble> EnvelopeProblemSOS<IPMDouble>::construct_SOS_instance() {
    unsigned const NUM_POLYNOMIALS = _polynomials_bounds.size();
    unsigned const VECTOR_LENGTH = _U;

    if (NUM_POLYNOMIALS == 0) {
        _logger->error("Please provide a polynomial.");
        exit(1);
    }
    if (NUM_POLYNOMIALS == 1) {
        _logger->warn("Instance trivial. Please detrivialize.");
        exit(1);
    }

    Constraints<IPMDouble> constraints;
    constraints.c = Vector::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH);
    constraints.c.block(0, 0, VECTOR_LENGTH, 1)
            = -_objectives_vector.template cast<IPMDouble>();

    constraints.A = Matrix::Zero((NUM_POLYNOMIALS - 1) * VECTOR_LENGTH, NUM_POLYNOMIALS * VECTOR_LENGTH);
    constraints.b = Vector::Zero((NUM_POLYNOMIALS - 1) * VECTOR_LENGTH);

    for (unsigned poly_idx = 0; poly_idx < NUM_POLYNOMIALS - 1; ++poly_idx) {
        Vector polynomial = _polynomials_bounds[poly_idx + 1].template cast<IPMDouble>();
        Matrix poly_block = Vector::Ones(VECTOR_LENGTH).asDiagonal();
        //corresponds to X variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, 0,
                            VECTOR_LENGTH, VECTOR_LENGTH)
                = -poly_block;
        //corresponds to Y_i variables
        constraints.A.block(poly_idx * VECTOR_LENGTH,
                            (poly_idx + 1) * VECTOR_LENGTH, VECTOR_LENGTH,
                            VECTOR_LENGTH)
                = poly_block;
        Vector casted_vec = _polynomials_bounds[0].template cast<IPMDouble>();
        constraints.b.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, 1)
                =
                polynomial - casted_vec;
    }

    _logger->info("Original SOS instance created.");
    if (_logger->level() == spdlog::level::trace) {
        constraints.print();
    }

    std::string num_prod_threads_str = getEnvVar("NUM_PRODUCT_THREADS");
    unsigned num_prod_threads = 1;
    if(num_prod_threads_str != ""){
        num_prod_threads = std::stoi(num_prod_threads_str);
    }

    //Construct Barrier function

    ProductBarrier<IPMDouble> *productBarrier = new ProductBarrier<IPMDouble>(num_prod_threads);

    for (unsigned poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
        auto sos_barrier = new InterpolantDualSOSBarrier<IPMDouble>(_d, _n);
        sos_barrier->configure(_config.get_child("InterpolantBarrier"));

        auto sum_barrier = new SumBarrier<IPMDouble>(_U);
        sum_barrier->add_barrier(sos_barrier);

        // Weighted barrier
        if (_use_weighted_polynomials) {
            Vector weighted_vec(3);
            //Add univariate polynomial 1 - x^2 which is non-negative for domain [-1,1]
            weighted_vec << 1., 0., -1.;
            auto sos_weighted_barrier = new InterpolantDualSOSBarrier<IPMDouble>(_d, weighted_vec);
            sos_weighted_barrier->configure(_config.get_child("InterpolantBarrier"));
            sum_barrier->add_barrier(sos_weighted_barrier);
        }

        productBarrier->add_barrier(sum_barrier);
    }

    auto dual_constraints = constraints.dual_system();

    Instance<IPMDouble> instance;
    instance.constraints = dual_constraints;
    instance.barrier = productBarrier;

    _logger->info("Dual formulation created.");
    if (_logger->level() == spdlog::level::trace) {
        instance.constraints.print();
    }
    return instance;
}


template <typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::print_solution(Solution<IPMDouble> sol) {
    // Note: the solution we are looking for is the RHS for the first polynomial - the SOS
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
}

template <typename IPMDouble>
void EnvelopeProblemSOS<IPMDouble>::plot_polynomials_and_solution(const Solution<IPMDouble> &sol) {

    std::cout << "Create picture of solution. Saved in plot.png..." << std::endl;
    int num_points = 1000;
    assert(num_points > 1);
    std::vector<double> x(num_points);
    assert(_hyperRectangle.size() == 1);

    IPMDouble x_min = _hyperRectangle[0].first;
    IPMDouble x_max = _hyperRectangle[0].second;

    IPM_DOUBLE const delta_x = x_max - x_min;
    x_min -= .05 * delta_x;
    x_max += .05 * delta_x;

    for (int j = 0; j < num_points; ++j) {
        IPMDouble d = x_min + j * (x_max - x_min) / (num_points - 1);
        Double dummy_D;
        x[j] = static_cast<Double>(d);
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

    if (_input_in_interpolant_basis) {
        //Have to be computed now, as we delayed computation in initialisation (was not necessary to create the instance)
        calculate_basis_polynomials();
    }

    InterpolantMatrix Q_interp = get_transformation_matrix();
    Matrix dummy_matrix;
//    Matrix Q = InterpolantMatrixToMatrix(Q_interp,dummy_matrix);

    std::vector<std::vector<double> > plots(poly_plots.size());
    for (unsigned poly_idx = 0; poly_idx < plots.size(); poly_idx++) {
        plots[poly_idx].resize(num_points);
        InterpolantVector poly_in_orig_basis(_U);
        poly_in_orig_basis = Q_interp * poly_plots[poly_idx];
        for (int i = 0; i < num_points; ++i) {
            //Evaluation of vector polynomial at a certain point. Extract method.
            assert(poly_in_orig_basis.size() > 0);
            InterpolantDouble eval = poly_in_orig_basis[0];
            for (int j = 1; j < poly_in_orig_basis.size(); ++j) {
                eval += poly_in_orig_basis[j] * pow(x[i], j);
            }
            plots[poly_idx][i] = static_cast<Double>(eval);
        }
    }

    const int w = 2000;
    const int h = w * 2 / 3;
    plt::figure_size(w, h);

    std::vector<double> &envelope_plot = plots[plots.size() - 1];
    std::vector<double> offset_envelope(envelope_plot.size());

    //Obtain plotted min and max.

    double y_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::min();

    for (unsigned i = 0; i < plots[0].size(); ++i) {
        if (x[i] < _hyperRectangle[0].first or x[i] > _hyperRectangle[0].second) {
            continue;
        }
        double local_y_min = std::numeric_limits<double>::max();
        for (unsigned plt_idx = 0; plt_idx < poly_plots.size() - 1; ++plt_idx) {
            y_min = std::min(y_min, plots[plt_idx][i]);
            local_y_min = std::min(local_y_min, plots[plt_idx][i]);
        }
        y_max = std::max(y_max, local_y_min);
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

    std::map<std::string, std::string> avx_keywords_string;
    avx_keywords_string["linestyle"] = "--";
    avx_keywords_string["color"] = "black";

    std::map<std::string, double> avx_keywords_double;
    avx_keywords_double["alpha"] = .5;
    plt::axvline(_hyperRectangle[0].first, 0, 1, avx_keywords_string, avx_keywords_double);
    plt::axvline(_hyperRectangle[0].second, 0, 1, avx_keywords_string, avx_keywords_double);

    // Plot a line whose name will show up as "log(x)" in the legend.
    // Add graph title
    std::string title_string = "Lower envelope";
    title_string += _use_weighted_polynomials ? ", weighted" : ", unweighted";
    title_string += ", degree " + std::to_string(_U - 1);
    title_string += ".";
    plt::title(title_string);
    // Enable legend.
    plt::legend();
    plt::save("plot");

    std::cout << "Done." << std::endl;
}

template <typename IPMDouble>
InterpolantMatrix EnvelopeProblemSOS<IPMDouble>::get_transformation_matrix() {
    InterpolantMatrix Q(_basis_polynomials.size(), _basis_polynomials.size());
    for (int i = 0; i < Q.rows(); ++i) {
        for (int j = 0; j < Q.rows(); ++j)
            Q(i, j) = _basis_polynomials[j][i];
    }
    return Q;
}
