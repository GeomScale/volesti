//
// Created by Bento Natura on 11/06/2020.
//

#include "EnvelopeProblemSDP.h"

EnvelopeProblemSDP::EnvelopeProblemSDP(int num_variables, int max_degree, HyperRectangle &hyperRectangle_) :
        _n(num_variables), _d(max_degree), monomialObject(num_variables, max_degree),
        _hyperRectangle(hyperRectangle_) {
    assert(num_variables == hyperRectangle_.size());
    construct_polynomial_matrix();
    construct_objective_matrix();
}

Matrix EnvelopeProblemSDP::get_objective_matrix() const {
    return _objectives_matrix;
}

void EnvelopeProblemSDP::add_polynomial(PolynomialSDP &polynomial) {
    _polynomials_bounds.push_back(polynomial);
}

PolynomialSDP EnvelopeProblemSDP::generate_zero_polynomial() {
    return Matrix::Zero(_objectives_matrix.rows(), _objectives_matrix.cols());
}

void EnvelopeProblemSDP::construct_polynomial_matrix() {
    monomialObject = MonomialsClass(_n, _d);

    auto monomials = monomialObject.getMonomials();
    _polynomial_product_matrix.resize(monomials.size());
    for (int i = 0; i < _polynomial_product_matrix.size(); ++i) {
        _polynomial_product_matrix[i].resize(monomials.size());
    }
    for (int i = 0; i < monomials.size(); ++i) {
        for (int j = 0; j < monomials.size(); ++j) {
            _polynomial_product_matrix[i][j] = MonomialsClass::prod(monomials[i], monomials[j]);
            monomialObject.print_human_readable(_polynomial_product_matrix[i][j]);
            std::cout << std::endl;
        }
    }
}
void EnvelopeProblemSDP::construct_objective_matrix() {
    _objectives_matrix.resize(_polynomial_product_matrix.size(), _polynomial_product_matrix.size());
    for (int i = 0; i < _polynomial_product_matrix.size(); ++i) {
        for (int j = 0; j < _polynomial_product_matrix.size(); ++j) {
            _objectives_matrix(i, j) = calculate_objective(_polynomial_product_matrix[i][j]);
        }
    }
    _objectives_matrix *= -1;
}

IPMDouble EnvelopeProblemSDP::calculate_objective(Monomial m) {
    return calculate_objective(m, 0);
}

IPMDouble EnvelopeProblemSDP::calculate_objective(Monomial m, int var) {
    assert(0 <= var and var <= _n);
    if (var == _n) {
        return 1;
    }
    IPMDouble exp = m[var] + 1;
    return 1. / exp * (pow(_hyperRectangle[var].second, exp) * calculate_objective(m, var + 1)
                       - pow(_hyperRectangle[var].first, exp) * calculate_objective(m, var + 1));
}


Instance EnvelopeProblemSDP::construct_SDP_instance() {
    unsigned const MATRIX_DIMENSION = _objectives_matrix.rows();
    unsigned const VECTOR_LENGTH = _objectives_matrix.rows() * _objectives_matrix.cols();
    unsigned const NUM_POLYNOMIALS = _polynomials_bounds.size();
    Constraints constraints;
    constraints.c = Vector::Zero((NUM_POLYNOMIALS + 1) * VECTOR_LENGTH);
    constraints.c.block(0, 0, VECTOR_LENGTH, 1) = MatrixToVector(_objectives_matrix);

    constraints.A = Matrix::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH, (NUM_POLYNOMIALS + 1) * VECTOR_LENGTH);
    constraints.b = Vector::Zero(NUM_POLYNOMIALS * VECTOR_LENGTH);

    for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
        PolynomialSDP polynomial = _polynomials_bounds[poly_idx];
        Matrix poly_block = Vector::Ones(VECTOR_LENGTH).asDiagonal();
        //corresponds to X variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, VECTOR_LENGTH) = poly_block;
        //corresponds to Y_i variables
        constraints.A.block(poly_idx * VECTOR_LENGTH, (poly_idx + 1) * VECTOR_LENGTH, VECTOR_LENGTH,
                            VECTOR_LENGTH) = poly_block;
        constraints.b.block(poly_idx * VECTOR_LENGTH, 0, VECTOR_LENGTH, 1) = MatrixToVector(polynomial);
    }

    //Construct Barrier function

    ProductBarrier *productBarrier = new ProductBarrier;
    auto X_barrier = new FullSpaceBarrier(VECTOR_LENGTH);
    productBarrier->add_barrier(X_barrier);

    for (int poly_idx = 0; poly_idx < NUM_POLYNOMIALS; ++poly_idx) {
        auto sdp_barrier = new SDPStandardBarrier(MATRIX_DIMENSION);
        productBarrier->add_barrier(sdp_barrier);
    }

    Instance instance;
    instance.constraints = constraints;
    instance.barrier = productBarrier;
    instance.constraints.print();

    return instance;
}

void EnvelopeProblemSDP::print_solution(const Solution &sol) {
    Vector v = sol.x.segment(0, _objectives_matrix.rows() * _objectives_matrix.cols());
    Matrix M = VectorToSquareMatrix(v, _objectives_matrix.rows());
    print_polynomial(M);
}

Monomial EnvelopeProblemSDP::get_matrix_entry(unsigned const row, unsigned const col) {
    assert(_polynomial_product_matrix.size() >= row);
    assert(_polynomial_product_matrix[row].size() >= col);
    return _polynomial_product_matrix[row][col];
}

void EnvelopeProblemSDP::print_polynomial(Matrix M) const {
    assert(M.cols() == _objectives_matrix.cols() and M.rows() == _objectives_matrix.rows());
    std::cout << "Polynomial is: " << std::endl;
    for (int row = 0; row < M.rows(); row++) {
        for (int col = 0; col <= row; col++) {
            IPMDouble coeff = (row == col) ? M(row, col) : M(row, col) + M(col, row);
            if (coeff) {
                std::cout << coeff << " * ";
                monomialObject.print_human_readable(_polynomial_product_matrix[row][col]);
                std::cout << std::endl;
            }
        }
    }
    std::cout << "--------------" << std::endl;
}

//Assume that monomial is ordered as 1, x, x*x, ...
IPMDouble EnvelopeProblemSDP::univariate_monomial_evaluation(Monomial const m, IPMDouble const x) {
    if (m[0] == 0) {
        return 1;
    }
    return pow(x, m[0]);
}

IPMDouble EnvelopeProblemSDP::univariate_polynomial_evaluation(PolynomialSDP const poly, IPMDouble x) {
    IPMDouble eval = 0.;
    for (int row = 0; row < poly.rows(); ++row) {
        for (int col = 0; col < poly.cols(); ++col) {
            eval += poly(row, col) * univariate_monomial_evaluation(_polynomial_product_matrix[row][col], x);
        }
    }
    return eval;
}

void EnvelopeProblemSDP::plot_polynomials_and_solution(const Solution &sol) {

    int num_points = 1000;
    assert(num_points > 1);
    std::vector<double> x(num_points);
    assert(_hyperRectangle.size() == 1);
    IPMDouble x_min = _hyperRectangle[0].first;
    IPMDouble x_max = _hyperRectangle[0].second;
    for (int j = 0; j < num_points; ++j) {
        IPMDouble d = x_min + j * (x_max - x_min) / (num_points - 1);
        x[j] = InterpolantDoubletoIPMDouble(d,x[j]);
    }

    std::vector<PolynomialSDP> poly_plots;
    for (PolynomialSDP &poly : _polynomials_bounds) {
        poly_plots.push_back(poly);
    }

    Vector v = sol.x.segment(0, _objectives_matrix.rows() * _objectives_matrix.cols());
    PolynomialSDP sol_poly = VectorToSquareMatrix(v, _objectives_matrix.rows());

    std::cout << "Objective of solution is: " << _objectives_matrix.cwiseProduct(sol_poly).sum() << std::endl;

    poly_plots.push_back(sol_poly);

    std::vector<std::vector<double> > plots(poly_plots.size());
    for (int poly_idx = 0; poly_idx < plots.size(); poly_idx++) {
        plots[poly_idx].resize(num_points);
        for (int i = 0; i < num_points; ++i) {
            //TODO: if plots fail, check here first
            IPMDouble ipm_d = univariate_polynomial_evaluation(poly_plots[poly_idx],x[i]);
            plots[poly_idx][i] = InterpolantDoubletoIPMDouble(ipm_d, plots[poly_idx][i]);
        }
    }

    plt::figure_size(1200, 780);
//        plt::xlim(x_min, x_max);
//        plt::ylim(-700, 500);

    for (int p_idx = 0; p_idx < poly_plots.size(); ++p_idx) {
        plt::plot(x, plots[p_idx]);
    }
    plt::named_plot("lower envelope", x, plots[plots.size() - 1]);

    // Plot a line whose name will show up as "log(x)" in the legend.
    // Add graph title
    plt::title("Lower envelope");
    // Enable legend.
    plt::legend();
    plt::save("plot");
}


