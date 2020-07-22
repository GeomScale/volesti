//
// Created by test Bento Natura on 22/07/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_INTERPOLANTDUALSOSBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_INTERPOLANTDUALSOSBARRIER_H

#include "LHSCB.h"

class InterpolantDualSOSBarrier : public LHSCB {

public:
    InterpolantDualSOSBarrier() : LHSCB() {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_) :
            InterpolantDualSOSBarrier(max_polynomial_degree_, Vector::Ones(1)) {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, Vector poly_g) : _max_polynomial_degree(
            max_polynomial_degree_) {

        assert(poly_g.rows() <= max_polynomial_degree_ + 1);

        _L = _max_polynomial_degree + 2 - poly_g.rows();
        _U = 2 * _max_polynomial_degree + 1;


        _preintermediate_matrix = Matrix(_U, _L);
        _intermediate_matrix = Matrix(_L, _L);
        _intermediate_LLT = Eigen::LLT<Matrix>(_L);
        _V = Matrix(_L, _U);
        _Q = Matrix(_V.cols(), _V.cols());

        _num_variables = _U;
        _unisolvent_basis.resize(_U);

        for (unsigned i = 0; i < _unisolvent_basis.size(); ++i) {
            BoostDouble cos_i = boost::multiprecision::cos(i * boost::math::constants::pi<BoostDouble>() / (_U - 1));
            InterpolantDouble dummy_ipm;
            InterpolantDouble cos_val = InterpolantDoubletoIPMDouble(cos_i, dummy_ipm);
            _unisolvent_basis[i] = cos_val;
        }

//        InterpolantMatrix P_interp = InterpolantMatrix(_U, _L);

        //TODO: Figure out how choice of P could influence condition / stability of maps.

        //Use monomial standard basis to orthogonalize

        //_P is used in the Interior Point Method. Therefore we need to convert the multi-precision
        // floating-point into the IPM floating point precision

        //TODO: This interpolant Matrix only needs to be found once and can then be reused;
        std::cout << "Construct interpolant point Matrix P..." << std::endl;

        //Construct P_interp via Chebyshev basis.

//        for (unsigned row = 0; row < P_interp.rows(); ++row) {
//            for (unsigned col = 0; col < P_interp.cols(); ++col) {
//                //TODO: make more efficient.
//                P_interp(row, col) = col == 0 ? InterpolantDouble(1.) : pow(_unisolvent_basis[row], col);
//            }
//        }

        //Alternative approach of finding _P via Chebyshev basis
        Eigen::MatrixXd cheb_P = ChebTools::u_matrix_library.get(_U - 1).block(0, 0, _U, _L);

        _g = Vector::Zero(_U);
        for (int p = 0; p < _U; ++p) {
            _g(p) = poly_g(0);
            for (int i = 1; i < poly_g.rows(); i++) {
                _g(p) += poly_g(i) * pow(_unisolvent_basis[p], i).convert_to<IPMDouble>();
            }
        }
        _g_g_transpose = _g * _g.transpose();

        //TODO: Option to compute cheb_P via InterpolantDouble;

        //TODO: Figure out whehter orthogonalization could be done in double precision to speed up initialisation.
        std::cout << "Constructed." << std::endl;
        std::cout << "Orthogonalize..." << std::endl;
        cxxtimer::Timer orth_timer;
        orth_timer.start();
        Matrix P_tmp = cheb_P.cast<IPMDouble>();
        Matrix P_ortho = P_tmp.householderQr().householderQ();
        P_ortho.colwise().hnormalized();
        _P = P_ortho.block(0, 0, _U, _L);

        orth_timer.stop();
        std::cout << "Orthogonalization done in " << orth_timer.count<std::chrono::milliseconds>() / 1000.
                  << " seconds." << std::endl;

    };

    bool update_gradient_hessian_LLT(Vector x, bool check_interior_only = false);

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Eigen::LLT<Matrix> llt(Vector x, bool symmetrize = 0) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    //TODO: better solution for implementation concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

    std::vector<InterpolantDouble> &get_basis() {
        return _unisolvent_basis;
    }

    Matrix get_P() {
        return _P;
    }

private:
    unsigned _max_polynomial_degree;
    std::vector<InterpolantDouble> _unisolvent_basis;
    Matrix _intermediate_matrix;
    Matrix _preintermediate_matrix;
    Eigen::LLT<Matrix> _intermediate_LLT;
    Matrix _Q;
    Matrix _V;
    unsigned _L, _U;

    //Weighted polynomials
    Vector _g;
    Matrix _g_g_transpose;
    Matrix _P;
};


#endif //NONSYMMETRICCONICOPTIMIZATION_INTERPOLANTDUALSOSBARRIER_H
