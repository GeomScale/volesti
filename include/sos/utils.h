// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_UTILS_H
#define SOS_UTILS_H

//Preprocessor directive allows us to forbid Eigen to allocate memory. Temporariliy helps to debug where allocation might slow down the program.
#define EIGEN_RUNTIME_NO_MALLOC

//Note MKL Macro is set in CmakeLists file.

//#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cxxtimer.hpp>
#include "spdlog/spdlog.h"
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "ChebTools/ChebTools.h"


#ifndef DIGITS_PRECISION
#define DIGITS_PRECISION 50
#endif

typedef boost::multiprecision::cpp_dec_float<DIGITS_PRECISION> mp_backend;
typedef boost::multiprecision::number<mp_backend, boost::multiprecision::et_on> BoostDouble;
//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<DIGITS_PRECISION> > BoostDouble;
typedef BoostDouble InterpolantDouble;

typedef long double long_double;

#ifdef IPM_DOUBLE
typedef IPM_DOUBLE Double;
#else
typedef double Double;
#endif

//Change typedef here to use different double type in interior point method.
#ifdef IPM_USE_DOUBLE
//typedef Double IPMDouble;
#else
typedef BoostDouble IPMDouble;
#endif

typedef Eigen::Matrix<BoostDouble, Eigen::Dynamic, Eigen::Dynamic> BoostMatrix;
typedef Eigen::Matrix<BoostDouble, Eigen::Dynamic, 1> BoostVector;

typedef Eigen::Matrix<InterpolantDouble, Eigen::Dynamic, Eigen::Dynamic> InterpolantMatrix;
typedef Eigen::Matrix<InterpolantDouble, Eigen::Dynamic, 1> InterpolantVector;

typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> DoubleMatrix;
typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> DoubleVector;

//Note: Boost Dependency
namespace pt = boost::property_tree;

template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

//Stack the columns of a square m x m  matrix to a vector of length m x m.
template<typename T>
Vector<T> StackMatrixToVector(Matrix<T> M) {
    assert(M.rows() == M.cols());
    Eigen::Map<Matrix<T>> x(M.data(), M.rows() * M.cols(), 1);
    return x;
}

//Unstack vector
template<typename T>
Matrix<T> UnstackVectorToMatrix(Vector<T> v, unsigned matrix_dimension) {
    assert(v.rows() == matrix_dimension * matrix_dimension);
    Eigen::Map<Matrix<T> > M(v.data(), matrix_dimension, matrix_dimension);
    return M;
}

//template <typename T>
//class Constraints<T>;

template<typename T>
class Solution {
public:
    template<typename U>
    Solution<U> cast() {
           Solution<U> sol;
           sol.x = x.template cast<U>();
           sol.s = s.template cast<U>();
           //Note: Boost Dependency.
           sol.centrality = boost::numeric_cast<U>(centrality);
           sol.gap= boost::numeric_cast<U>(gap);
           return sol;
    }

    Vector<T> x;
    Vector<T> s;
    T centrality;
    T gap;
};

//TODO: Need full row rank matrices for IPM. Also, is preprocessing A, e.g. row-echelon form useful?
template<typename T>
class Constraints {
public:
    Matrix<T> A;
    Vector<T> b;
    Vector<T> c;

    Constraints() {};

    Constraints(Matrix<T> A_, Vector<T> b_, Vector<T> c_) : A(A_), b(b_), c(c_) {};

    void print() {
        std::cout << "Constraints are as follows. Constraint matrix is A: " << std::endl;
        std::cout << A << std::endl;
        std::cout << "Objective c is " << std::endl;
        std::cout << c.transpose() << std::endl;
        std::cout << "RHS b is " << std::endl;
        std::cout << b.transpose() << std::endl;
    }

    //TODO: use optional argument to indicate sparseness.
    Constraints dual_system() {
        std::cout << "Create dual system... " << std::endl;
        Constraints dual_constraints;
//        dual_constraints.c = A.colPivHouseholderQr().solve(b);

        //TODO: use proper tolerance / reference.

        Eigen::SparseMatrix<T> A_top_sparse = A.transpose().sparseView(0, 1e-10);
        Eigen::SparseMatrix<T> A_sparse = A.sparseView(0, 1e-10);

        A_top_sparse.makeCompressed();
        A_sparse.makeCompressed();

        Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int> > QR_top_sparse;
        Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int> > QR_sparse;

        QR_top_sparse.compute(A_top_sparse);
        QR_sparse.compute(A_sparse);

        dual_constraints.c = QR_sparse.solve(b);

        Matrix<T> QR_from_sparse(QR_top_sparse.matrixQ());
//        Matrix QR = A.transpose().householderQr().householderQ();

        dual_constraints.A = QR_from_sparse.block(0, A.rows(), QR_from_sparse.rows(),
                                                  QR_from_sparse.cols() - A.rows()).transpose();
        dual_constraints.b = dual_constraints.A * c;

        std::cout << "Done." << std::endl;

        //TODO: use different measure to calculate centrality error
        assert((dual_constraints.A * A.transpose()).norm() < 10e-5);
        return dual_constraints;
    }
};

//Naive implementation
class DegreeTuple {
public:
    DegreeTuple(const int num_vars, const unsigned max_degree_) {
        max_degree = max_degree_;
        v.resize(num_vars, 0);
    }

    bool next() {
        for (int i = v.size() - 1; i >= 0; i--) {
            if (v[i] < max_degree) {
                v[i]++;
                return true;
            }
            v[i] = 0;
        }
        return false;
    }

    bool valid() {
        int sum = 0;
        for (auto el : v) {
            sum += el;
        }
        return sum <= max_degree;
    }

    bool next_valid() {
        if (!next()) {
            return false;
        }
        while (not valid()) {
            if (not next()) {
                return false;
            }
        }
        return true;
    }

    std::vector<unsigned> &get_tuple() {
        return v;
    }

    void print_tuple() {
        for (auto el : v) {
            std::cout << el << " ";
        }
        std::cout << std::endl;
    }

private:
    unsigned max_degree;
    std::vector<unsigned> v;
};

class AllCombinationTuple {
public:
    AllCombinationTuple(std::vector<unsigned> const bounds_) {

        bounds = bounds_;
        v.resize(bounds.size(), 0);
    }

    bool next() {
        for (int i = v.size() - 1; i >= 0; i--) {
            if (bounds[i] > v[i]) {
                v[i]++;
                return true;
            }
            v[i] = 0;
        }
        return false;
    }

    std::vector<unsigned> &get_combination() {
        return v;
    }

    std::vector<unsigned> bounds;
    std::vector<unsigned> v;
};


//Implementation copied from https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

template<typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}

//Access environment variables

std::string getEnvVar(std::string const &key);

template<typename T>
class OrthogonaPMatrixLibrary {
private:
    std::map<std::pair<int, int>, Matrix<T> > matrices;

    void build(int L, int U) {
        Eigen::MatrixXd cheb_P = ChebTools::u_matrix_library.get(U - 1).block(0, 0, U, L);
        Matrix<T> P_tmp = cheb_P.cast<T>();
        Matrix<T> P_ortho = P_tmp.householderQr().householderQ();
        P_ortho.colwise().hnormalized();
        matrices[std::pair<int, int>(L, U)] = P_ortho.block(0, 0, U, L);
    }

public:
    /// Get the \f$\mathbf{U}\f$ matrix of degree N
    const Matrix<T> &get(int L, int U) {
        auto it = matrices.find(std::pair<int, int>(L, U));
        if (it != matrices.end()) {
            return it->second;
        } else {
            build(L, U);
            return matrices.find(std::pair<int, int>(L, U))->second;
        }
    }
};

template<typename T>
static OrthogonaPMatrixLibrary<T> orthogonal_P_Matrix_library;


template<typename _MatrixType, int _UpLo> class CustomLLT : public Eigen::LLT<_MatrixType, _UpLo> {
public:

    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;

    CustomLLT<_MatrixType, _UpLo>() :  Eigen::LLT<_MatrixType, _UpLo>() {};

    CustomLLT<_MatrixType, _UpLo>(unsigned int n) :  Eigen::LLT<_MatrixType, _UpLo>(n) {};

    void copy_and_scale(const CustomLLT<_MatrixType, _UpLo> & other, Scalar scalar){
        this->m_matrix = other.m_matrix * scalar;
        this->m_l1_norm = other.m_l1_norm;
        this->m_isInitialized = other.m_isInitialized;
        this->m_info = other.m_info;
    }
};

#endif //SOS_UTILS_H



