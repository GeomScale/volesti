//
// Created by panagiotis on 25/6/2019.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

//#include "Eigen"
#include "../../Spectra/SymGEigsSolver.h"
#include "../../Spectra/MatOp/DenseSymMatProd.h"
#include "../../Spectra/MatOp/SparseCholesky.h"
#include <vector>
//#include <Eigen/Eigen>
#include <chrono>
#include <limits>
#include "optimization/SpectraLU.h"
#include "SpectraLU.h"

//const double ZERO = 0.0000000001;
//double ORACLE_TIME;
//double REFLECTION_TIME;
//int BOUNDARY_CALLS;

//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MT_ROWMAJOR;
//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
//typedef Eigen::SparseMatrix<double> SpMat;

/**
 * A linear matrix inequality A_0 + sum(x * F_i)
 */
template <class TMT, class TVT>
class LMI {

    TMT A0;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MT_ROWMAJOR2;
    // the matrices A_i, i>0
    std::vector<TMT> matrices;
    MT_ROWMAJOR2 vectorMatrix;
    MT_ROWMAJOR2 gradientMatrix;
    TVT a;
    typedef typename std::vector<TMT>::iterator Iter;

public:

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MT_ROWMAJOR;
    typedef TMT MT;
    typedef TVT VT;

    LMI() {};

    LMI(std::vector<MT> &matrices) {
        this->A0 = matrices[0];
        for (int i = 1; i < matrices.size(); i++) {
            this->matrices.push_back(matrices[i]);
        }

        setGradientMatrix();
        setVectorMatrix();
    }

    LMI(const LMI &lmi) {
        this->A0 = lmi.A0;
        this->matrices = lmi.matrices;
        setVectorMatrix();
        setGradientMatrix();
    }

    void set_Ai(const MT &Ai, const int i) {
        matrices[i] = Ai;
    }

    void set_A0(const MT &MatA0) {
        A0 = MatA0;
    }

    /**
     * Evaluate the lmi for vector x
     *
     * @param x
     * @return
     */
    void evaluate(const VT &x, MT &res) {
        int dim = A0.rows();
        res = MT::Zero(dim, dim);
        res += A0;
        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res.noalias() += x(i) * (*iter);
    }

    void evaluate_revised(const VT &x, MT &res) {
        evaluateWithoutA0_revised(x, res);
        res += A0;
    }


    bool isSingular(VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() == 0)
                return true;

        return false;
    }

    bool isSingular(VT &x, double approx) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (fabs(eivals(i).real()) <= fabs(approx))
                return true;

        return false;
    }

    bool isNegativeDefinite(const VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() >= 0)
                return false;

        return true;
    }

    bool isNegativeSemidefinite(const VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename  Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++) {
//            std::cout << eivals(i).real() << "\n";
            if (eivals(i).real() > 0)
                return false;
        }

        return true;
    }

    bool isPositiveSemidefinite(VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() < 0)
                return false;

        return true;
    }

    bool isPositiveDefinite(VT &x) {
        MT mt;
        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() <= 0)
                return false;

        return true;
    }


    bool isPositiveSemidefinite(VT &x, MT &mt, VT &minEigenVector, double &minEigenvalue) {
        bool res = true;

        evaluate(x, mt);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
        minEigenVector.resize(solver.eigenvectors().col(0).rows());

        for (int i = 0; i < minEigenVector.rows(); i++) {
            minEigenVector(i) = solver.eigenvectors().col(0)(i).real();

        }

        double min = eivals(0).real() / minEigenVector.norm();
        int minIndex = 0;

        for (int i = 0; i < eivals.rows(); i++) {
            for (int j = 0; j < minEigenVector.rows(); j++) {
                minEigenVector(j) = solver.eigenvectors().col(i)(j).real();
            }

            if (eivals(i).real() < 0)
                res = false;

            if (eivals(i).real() / minEigenVector.norm() < min) {
                min = eivals(i).real() / minEigenVector.norm();
                minIndex = i;
            }
        }


        minEigenvalue = eivals(minIndex).real() / minEigenVector.norm();

        for (int i = 0; i < minEigenVector.rows(); i++) {
            minEigenVector(i) = solver.eigenvectors().col(minIndex)(i).real();

        }
        minEigenVector.normalize();

        return res;
    }

    void evaluateHMC(const VT &x, MT &res) const {
        int dim = A0.rows();
        res = MT(A0);
        int i = 0;

        for (auto M : matrices)
            res += x(i++) * M;
    }

    int getDim() const {
        return matrices.size();
    }

    int getMatricesDim() const {
        return A0.rows();
    }

    const std::vector<MT> &getMatrices() const {
        return matrices;
    }

    /**
     * Evaluate the lmi for vector x without taking int account matrix A0
     *
     * @param x
     * @return
     */
    void evaluateWithoutA0(const VT &x, MT &res) {
        long dim = A0.rows();

        res = MT::Zero(dim, dim);

        int i = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++)
            res.noalias() += x(i) * (*iter);
    }

    /**
 * Evaluate the lmi for vector x without taking int account matrix A0
 *
 * @param x
 * @return
 */
    void evaluateWithoutA0_revised(const VT &x, MT &res) {

        long dim = A0.rows();

        a = vectorMatrix * x;

        double *data = res.data();
        double *v = a.data();

        int at = 0;

        // copy lower triangular
        for (int at_col = 0; at_col < dim; at_col++) {
            int col_offset = at_col * dim;
            double *target = data + col_offset + at_col;

            for (int at_row = at_col; at_row < dim; at_row++) {
                *(target++) = *(v++);
            }
        }

        v = a.data();

        // copy upper triangular
        for (int at_row = 0; at_row < dim; at_row++) {
            double *target = data + at_row + at_row * dim;

            for (int at_col = at_row; at_col < dim; at_col++) {
                *target = *(v++);
                target = target + dim;
            }
        }

    }

    void setVectorMatrix() {
        int dim = getMatricesDim();
        int newDim = dim * (dim + 1) / 2;

        a.resize(newDim);

        vectorMatrix.resize(newDim, getDim());

        int atMatrix = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, atMatrix++) {
            int i = 0;

            for (int at_row = 0; at_row < dim; at_row++)
                for (int at_col = at_row; at_col < dim; at_col++) {
                    vectorMatrix(i++, atMatrix) = (*iter)(at_row, at_col);
                }

        }

    }

    void setGradientMatrix() {
        int m = getMatricesDim();
        int d = getDim();

        gradientMatrix.resize(d*m, m);
        int atMatrix = 0;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, atMatrix++) {
            gradientMatrix.block(atMatrix*m,0,m,m) = matrices[atMatrix];

        }

    }

    const MT_ROWMAJOR &getGradientMatrix() const {
        return gradientMatrix;
    }

    const MT &getAi(int i) const {
        return matrices[i];
    }

    const std::vector<MT> getMatrices(){
        return matrices;
    }

    const MT &getA0() const {
        return A0;
    }

    void setA0(MT &A0) {
        this->A0 = A0;
    }

    void addMatrix(MT &matrix) {
        matrices.push_back(matrix);
    }

    void print() {
        std::cout << "F0" << "\n" << A0 << "\n";
        int i = 1;

        for (Iter iter = matrices.begin(); iter != matrices.end(); iter++, i++) {
            std::cout << "F" << i << "\n";
            std::cout << *iter << "\n";
        }
    }

    void changeInequalityDirection() {
        A0 = -1 * A0;

        for (int i = 0; i < matrices.size(); i++)
            matrices[i] = -1 * matrices[i];
    }

};

template <class TLMI, class Point>
class Spectrahedron {
    /**
     * The collection of matrices that constitute the linear matrix
     * inequality describing the spectrahedron
     */
    

    TLMI lmi;
    typename TLMI::VT U;//for computing gradient
    typename TLMI::MT LMIatP;
    double maxDouble = std::numeric_limits<double>::max();
    double minDouble = std::numeric_limits<double>::lowest();

public:

    typedef Point PolytopePoint;
    typedef TLMI LMI;
    typedef typename TLMI::MT MT;
    typedef typename TLMI::VT VT;
    typedef typename TLMI::MT_ROWMAJOR MT_ROWMAJOR;
    
    Spectrahedron() {}

    Spectrahedron(const Spectrahedron &spectrahedron) {
        LMI lmi;
        this->lmi = LMI(spectrahedron.getLMI());
        U.resize(lmi.getDim()*lmi.getMatricesDim());
        LMIatP.resize(getLMI().getMatricesDim(), getLMI().getMatricesDim());
    }

    Spectrahedron(LMI &lmi) {
        this->lmi = lmi;
        U.resize(lmi.getDim()*lmi.getMatricesDim());
        LMIatP.resize(getLMI().getMatricesDim(), getLMI().getMatricesDim());
    }

    const LMI &getLMI() const {
        return lmi;
    }

    int dimension() {
        return getLMI().getDim();
    }

    void print_matrices() {
        std::cout<<"A0 = \n"<<getLMI().getA0()<<"\n"<<std::endl;
        for (int i = 0; i < dimension(); ++i) {
            std::cout<<"A"<<i+1<<" = \n"<<getLMI().getAi(i)<<"\n"<<std::endl;
        }
    }

    /**
     * Compute the intersection of a 1D line and the spectrahedron by finding the
     * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
     *
     * @param position
     * @param direction
     * @return (minimum positive eigenvalue, maximum negative eigenvalue)
     */
    std::pair<double, double> boundaryOracle(const VT &position, const VT &direction) {
        MT A;
        lmi.evaluate(position, A);
        MT B;
        lmi.evaluateWithoutA0(-1 * direction, B);


        Eigen::GeneralizedEigenSolver<MT> ges(A, B);
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;


        for (int i = 0; i < alphas.rows(); i++) {
            if (betas(i) == 0) //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }

        // for numerical stability
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive == maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
        if (lambdaMaxNegative == minDouble) lambdaMaxNegative = 0;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /**
    * Compute the intersection of a 1D line and the spectrahedron by finding the
    * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
    *
    * Take also in account the halfspace ax<=b
     *
    * @param position
    * @param direction
    * @return (minimum positive eigenvalue, maximum negative eigenvalue)
    */
    std::pair<double, double>
    boundaryOracleEfficient(const VT &position, const VT &direction, const VT &a, double b) {
        MT A;
        lmi.evaluate(position, A);
        MT B;
        lmi.evaluateWithoutA0(-1 * direction, B);
        //BOUNDARY_CALLS++;
        // Construct matrix operation object using the wrapper classes
        Spectra::DenseSymMatProd<double> op(B);
        Spectra::DenseCholesky<double> Bop(-A);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 3);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();
        // Retrieve results
        Eigen::VectorXd evalues;
//        Eigen::MatrixXd evecs;
        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
//            evecs = geigs.eigenvectors();
        }

        double lambdaMaxNegative;
        double lambdaMinPositive;

        if (nconv == 2) {
            lambdaMaxNegative = -1 / evalues(0);
            lambdaMinPositive = -1 / evalues(1);
        } else {
            lambdaMaxNegative = lambdaMinPositive = 0;
        }

        // check the cutting plane
        double lambda = (b - a.dot(position)) / a.dot(direction);
        if (lambda > 0 && lambda < lambdaMinPositive)
            lambdaMinPositive = lambda;
        if (lambda < 0 && lambda > lambdaMaxNegative)
            lambdaMaxNegative = lambda;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /**
    * Compute the intersection of a 1D line and the spectrahedron by finding the
    * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
    *
    * Take also in account the halfspace ax<=b
     *
    * @param position
    * @param direction
    * @return (minimum positive eigenvalue, maximum negative eigenvalue)
    */
    std::pair<double, double> boundaryOracle(const VT &position, const VT &direction, const VT &a, double b) {
        MT A;LMIatP.resize(getLMI().getMatricesDim(), getLMI().getMatricesDim());
        lmi.evaluate(position, A);
        MT B;
        lmi.evaluateWithoutA0(-1 * direction, B);

        Eigen::GeneralizedEigenSolver<MT> ges(A, B);
        //BOUNDARY_CALLS++;
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }


        // for numerical stability
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive == maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
        if (lambdaMaxNegative == minDouble) lambdaMaxNegative = 0;


        // check the cutting plane
        double lambda = (b - a.dot(position)) / a.dot(direction);
        if (lambda > 0 && lambda < lambdaMinPositive)
            lambdaMinPositive = lambda;
        if (lambda < 0 && lambda > lambdaMaxNegative)
            lambdaMaxNegative = lambda;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    void print() {
        this->lmi.print();
    }

    void changeInequalityDirection() {
        lmi.changeInequalityDirection();
    }

    bool isSingular(VT &x) {
        return lmi.isSingular(x);
    }


    bool isSingular(VT &x, double approx) {
        return lmi.isSingular(x, approx);
    }

    void set_Ai(const MT &Ai, const int i) {
        lmi.set_Ai(Ai, i);
    }

    //int dimension() {
        //return lmi.getDim();
    //}

    // Cholesky decomposition
    // It is important to be optimized (oxi axreiastoi upologismoi)
    //template <class MT>
    int is_in(Point &p) {

        lmi.evaluate_revised(p.getCoefficients(), LMIatP);
        //LMIatP += getLMI().getA0();
        //Check and put your code here [for Panagiotis)
        Spectra::DenseCholesky<double> Bop(-LMIatP);
        if (Bop.info() == Spectra::SUCCESSFUL) {
            return -1;
        }
        return 0;
    }

    //template<class Point>
    class BoundaryOracleCDHRSettings {
    public:
        MT LMIatP;
        double max_segment_length;
        bool first; //true if first call of the boundary oracle
        bool hasCuttingPlane;

        BoundaryOracleCDHRSettings(const int& dim) {
            LMIatP.resize(dim, dim);
            first = true;
            max_segment_length = 0;
            hasCuttingPlane = false;
        }

        void setMaxSegmentLength(double minLambdaPositive, double maxLambdaNegative) {
            double lambda = minLambdaPositive - maxLambdaNegative;

            if (lambda > max_segment_length)
                max_segment_length = lambda;
        }

        void resetMaxSegmentLength() {
            max_segment_length = 0;
        }

        double maxSegmentLength() {
            return max_segment_length;
        }
    };

    //template<class Point>
    std::pair<double, double>
    boundaryOracleCDHR(const VT &position, const int& coordinate, const VT &a, const double &b,
                           BoundaryOracleCDHRSettings &settings) {

        //ORACLE_TIME++;
        if (settings.first) {
            lmi.evaluate_revised(position, settings.LMIatP);
        }

        // Construct matrix operation object using the wrapper classes

        Spectra::DenseSymMatProd<double> op(lmi.getMatrices()[coordinate]);
        Spectra::DenseCholesky<double> Bop(-settings.LMIatP);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 15 < lmi.getMatricesDim() ? 15 : lmi.getMatricesDim());

        // Initialize and compute
        geigs.init();
        //std::cout<<"CDHR initialization done"<<std::endl;
        int nconv = geigs.compute();
        //std::cout<<"CDHR computation done"<<std::endl;


        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::MatrixXd evecs;

        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
            evecs = geigs.eigenvectors();
        }

        double lambdaMinPositive, lambdaMaxNegative;

        if (nconv == 2) {
            lambdaMinPositive = 1 / evalues(0);
            lambdaMaxNegative = 1 / evalues(1);
        } else {
            return {0,0};
        }

        if (lambdaMinPositive < 0)
            throw("Unbounded");

        if (lambdaMaxNegative > 0 || (lambdaMaxNegative == 0 && lambdaMinPositive ==0))
            throw("Unbounded");

        if (settings.hasCuttingPlane) {
            //std::cout<<"check cutting plane"<<std::endl;
            // check the cutting plane
            double lambda = (b - a.dot(position)) / a[coordinate];

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }

        settings.setMaxSegmentLength(lambdaMinPositive, lambdaMaxNegative);

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    //template<class Point>
    class BoundaryOracleBilliardSettings {
    public:
        MT B, LMIatP;
        Point genEigenvector;
        double max_segment_length;
        bool first; //true if first call of the boundary oracle
        bool computeB;
        MT Obj;
        bool isObjComputed;

        BoundaryOracleBilliardSettings(int dim) {
            first = true;
            LMIatP.resize(dim, dim);
            B.resize(dim, dim);
            max_segment_length = 0;
            genEigenvector = Point(dim);
            Obj.resize(dim, dim);
            isObjComputed = false;
            computeB = true;
        }

        void setMaxSegmentLength(double lambda) {
            if (lambda > 0 && lambda > max_segment_length)
                max_segment_length = lambda;
        }

        void resetMaxSegmentLength() {
            max_segment_length = 0;
        }

        double maxSegmentLength() {
            return max_segment_length;
        }
    };

    //template<class Point>
    std::pair<double, bool>
    boundaryOracleBilliard(const VT &position, const VT &direction, const VT &a, const double &b,
                           BoundaryOracleBilliardSettings &settings, bool check_cutting_plane = true) {

        if (settings.first) {
            lmi.evaluate_revised(position, settings.LMIatP);
        }

        if (!settings.isObjComputed) {
            settings.isObjComputed = true;
            lmi.evaluateWithoutA0_revised(a, settings.Obj);
        }

        //BOUNDARY_CALLS++;

        if (settings.computeB)
            lmi.evaluateWithoutA0_revised(direction, settings.B);

        settings.computeB = true;

        // check the cutting plane
        double lambdaMinPositive = 0;
        bool hitCuttingPlane = false;

        if (check_cutting_plane) {
            double lambda = (b - a.dot(position)) / a.dot(direction);

            if (lambda > 0) {
                lambdaMinPositive = lambda;
                hitCuttingPlane = true;

                Spectra::DenseCholesky<double> Bop(-settings.LMIatP - lambda * settings.B);

                if (Bop.info() == Spectra::SUCCESSFUL) {
                   // REFLECTION_TIME++;
                    settings.setMaxSegmentLength(lambdaMinPositive);
                    return {lambdaMinPositive, hitCuttingPlane};
                }

            }
        }



        // Construct matrix operation object using the wrapper classes
        Spectra::DenseSymMatProd<double> op(settings.B);
        Spectra::DenseCholesky<double> Bop(-settings.LMIatP);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        Spectra::SymGEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 1, 15 < lmi.getMatricesDim() ? 15 : lmi.getMatricesDim());

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::MatrixXd evecs;

        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
            evecs = geigs.eigenvectors();
        }


        if (nconv == 1) {
            double lambda = 1 / evalues(0);
            if (lambda < 0)
                throw("Unbounded");

            if (lambdaMinPositive == 0 || lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            else
                hitCuttingPlane = true;

            settings.genEigenvector.set(evecs.col(0));
        } else {
            lambdaMinPositive = 0;
        }

        // check the cutting plane
//        bool hitCuttingPlane = false;
//
//        if (check_cutting_plane) {
//            double lambda = (b - a.dot(position)) / a.dot(direction);
//
//            if (lambda > 0 && lambda < lambdaMinPositive) {
//                lambdaMinPositive = lambda;
//                hitCuttingPlane = true;
//            }
//        }

        settings.setMaxSegmentLength(lambdaMinPositive);


        return {lambdaMinPositive, hitCuttingPlane};
    }


    //template<class Point>
    class BoundaryOracleBoltzmannHMCSettings {
    public:
        MT B0, B1, B2;
        Point genEigenvector;
        Point gradient;
        bool first; //true if first call of the boundary oracle
        double epsilon; //when a point is this close to the boundary, consider it a boundary point
        bool isBoundaryPoint;

        BoundaryOracleBoltzmannHMCSettings() {
            first = true;
            epsilon = 0.0001;
        }
    };

    template<class SettingClass>
    double boundaryOracle_Boltzmann_HMC(const Point &_position, const Point &_direction, const Point &_objectiveFunction,
                                        const double &temp, SettingClass &settings) {

        const VT &position = _position.getCoefficients();
        const VT &direction = _direction.getCoefficients();
        const VT &objectiveFunction = _objectiveFunction.getCoefficients();
        //BOUNDARY_CALLS++;
        unsigned int matrixDim = lmi.getMatricesDim();
//        if (!lmi.isNegativeSemidefinite(position)) throw "out\n";
//            std::cout << objectiveFunction << "\n";
//        if (first) {
        lmi.evaluateHMC(position, settings.B0);
        lmi.evaluateWithoutA0(objectiveFunction, settings.B2);
//        }

        lmi.evaluateWithoutA0(direction, settings.B1);
        MT B2temp = settings.B2 / (-2 * temp);

        // create pencil matrix
        MT AA(2 *matrixDim, 2 * matrixDim);
        MT BB(2 * matrixDim, 2 * matrixDim);

        BB.block(matrixDim, matrixDim, matrixDim, matrixDim) = -1 * settings.B0;
        BB.block(0, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        BB.block(matrixDim, 0, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);
        BB.block(0, 0, matrixDim, matrixDim) = B2temp;

        AA.block(0, matrixDim, matrixDim, matrixDim) = settings.B0;
        AA.block(0, 0, matrixDim, matrixDim) = settings.B1;
        AA.block(matrixDim, 0, matrixDim, matrixDim) = settings.B0;
        AA.block(matrixDim, matrixDim, matrixDim, matrixDim) = MT::Zero(matrixDim, matrixDim);


//        double frobeniusNorm = 0;
//        for (int i=0 ; i<2*matrixDim ; i++)
//            for (int j=0 ; j<2*matrixDim ; j++)
//                frobeniusNorm += BB(i,j) * BB(i,j);
//        frobeniusNorm = std::sqrt(frobeniusNorm);

//        SpMat A =AA.sparseView();
//        SpMat B=BB.sparseView();
//    std::cout << AA << "\n\n\n";
//        std::cout << BB << "\n";

        // Construct matrix operation object using the wrapper classes
//        Spectra::SparseSymMatProd<double> op(A);
//        Spectra::DenseSymMatProd<double> op(BB);
//        Spectra::SparseRegularInverse<double> Bop(B);

//        Spectra::SpectraLU<double> Bop(AA);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
//        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEIGS_CHOLESKY>
//                geigs(&op, &Bop, 2, 4);

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
//        Spectra::SymGEigsSolver<double, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<double>, Spectra::SpectraLU<double>, Spectra::GEIGS_REGULAR_INVERSE>
//                geigs(&op, &Bop, 1, 3);

        // Initialize and compute
//        geigs.init();
//        int nconv = geigs.compute();
//         Retrieve results
//        Eigen::VectorXd evalues;
//        Eigen::MatrixXd evecs;
//        if(geigs.info() == Spectra::SUCCESSFUL)
//        {
//            evalues = geigs.eigenvalues();
//            evecs = geigs.eigenvectors();
//        }
//
        double lambdaMinPositive = 0;

//        if (nconv == 1) {
//            lambdaMinPositive = 1/evalues(0);
//            settings.genEigenvector = Point(evecs.col(0).segment(matrixDim, matrixDim));
//        } else {
//            lambdaMinPositive = 0;
//        }

//        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(AA, -1*BB);
//
//        bool _first = true;
//        int index;
//        const VT& eivals = es.eigenvalues();
//
//        for (int i=0 ; i<eivals.rows() ; i++) {
//            if (eivals(i) > 0 && (_first || eivals(i) < lambdaMinPositive)) {
//                lambdaMinPositive = eivals(i);
//                index = i;
//                _first = false;
//            }
//        }
//
//        MT vecs = es.eigenvectors();
//        genEivector = vecs.col(index).segment(matrixDim, matrixDim);





//        this->getLMI().print();
//        std::cout << position << "\n" << direction << "\n" << temp<<"\n";
        Eigen::GeneralizedEigenSolver<MT> ges(AA, -BB);

        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        lambdaMinPositive = maxDouble;
        int index = 0;

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0 || alphas(i).imag() != 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);
//            std::cout << lambda << " e\n";
            if (lambda > 0 && lambda < lambdaMinPositive) {
                lambdaMinPositive = lambda;
                index = i;
            }
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative = lambda;
        }
//        std::cout << objectiveFunction <<"\n";
//        std::cout << lambdaMinPositive << "eval\n";
//        std::cout << AA <<"\n";
//        std::cout << BB;

        typename Eigen::GeneralizedEigenSolver<MT>::EigenvectorsType eivecs = ges.eigenvectors();
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivec = eivecs.col(index);

        settings.genEigenvector = Point(matrixDim);

        for (int i = 0; i < matrixDim; i++)
            settings.genEigenvector.set_coord(i, eivec(matrixDim + i).real());




//        std::cout << eivecs.col(index)<<"evec\n";
        return lambdaMinPositive;
    }

    template<class SettingClass>
    void compute_reflection(SettingClass &settings, Point &direction) {
        std::vector<MT> matrices = lmi.getMatrices();
        int dim = matrices.size();
        settings.gradient = Point(dim);

        Point q(dim),w(dim);
        for (int i = 0; i < dim; i++) {
          q = settings.genEigenvector;
          w = Point(matrices[i] * q.getCoefficients());
          settings.gradient.set_coord(i, q.dot(w));
            //settings.gradient.set_coord(i, settings.genEigenvector.dot(
                    //(settings.genEigenvector.matrix_left_product(matrices[i]))));
        }

        settings.gradient.normalize();
        settings.gradient = -1 * settings.gradient;

        Point t = ((-2.0 * direction.dot(settings.gradient)) * settings.gradient);
        direction = t + direction;
        settings.gradient = -1 * settings.gradient;
    }

    /*void compute_reflection(const VT &genEivector, VT &direction, MT &C) {
        VT gradient;
        std::vector<MT> matrices = lmi.getMatrices();
        int dim = matrices.size();
        gradient = VT::Zero(dim);

        for (int i = 0; i < dim; i++) {
            gradient(i) = genEivector.dot((matrices[i]) * genEivector);
        }

//        Eigen::SelfAdjointEigenSolver<MT> es;
//        es.compute(C);
//        double product = 1;
//        VT evecs = es.eigenvalues();
//        for (int i=0 ; i<evecs.rows() ; i++) {
//            if (evecs(i)!= 0)
//                product *= evecs(i);
//        }

//        Eigen::FullPivLU<MT> lu_decomp(C);
//        auto rank = lu_decomp.rank();

//        std::cout << rank << "\n";
//        gradient *= product;

        gradient = -gradient;
        gradient.normalize();

        gradient = ((-2.0 * direction.dot(gradient)) * gradient);
        direction += gradient;
    }*/

    //template<class Point>
    void compute_reflection(const Point &genEivector, Point &direction, const Point &p) {

//        auto t1 = std::chrono::steady_clock::now();

//        const std::vector<MT> &matrices = lmi.getMatrices();
        int d = getLMI().getDim();
        int m = lmi.getMatricesDim();
        const VT& coeffs = genEivector.getCoefficients();

        U.noalias() = lmi.getGradientMatrix() * coeffs;
        Point gradient(d);

        for (int i = 0; i < d; i++) {
//            gradient.set_coord(i, genEivector.dot(genEivector.matrix_left_product(matrices[i])));
            gradient.set_coord(i, U.segment(i*m,m).dot(coeffs));
        }

//        gradient = -1 * gradient;
        gradient.normalize();

        gradient *= (-2.0 * direction.dot(gradient));
        direction += gradient;


    }

    template <typename NT>
    void ComputeInnerBall(NT &diam, NT &radius) {

        BoundaryOracleCDHRSettings CDHRsettings(getLMI().getMatricesDim());
        CDHRsettings.LMIatP = getLMI().getA0();
        diam = 0.0;
        radius = maxDouble;

        VT center(dimension());
        Point v(dimension());// v(dimension());
        NT bb=0.0;

        for (unsigned int i = 0; i < dimension(); ++i) {
            //std::cout<<"i = "<<i<<std::endl;
            v = Point(dimension());
            //v = VT::Zero(getLMI().getMatricesDim());
            //v(i) = 1.0;
            //std::pair<NT, NT> min_max = boundaryOracleCDHR(center, i, center, bb, CDHRsettings);
            v.set_coord(i, 1.0);
            std::pair<NT, NT> min_max = boundaryOracle(center, v.get_coefficients());
            //std::cout<<"min_max.first = "<<min_max.first<<", min_max.second = "<<min_max.second<<std::endl;
            if (min_max.first < radius) radius = min_max.first;
            if (-min_max.second < radius) radius = -min_max.second;
            if (min_max.first-min_max.second > diam ) diam = min_max.first-min_max.second;
        }

        radius = radius / std::sqrt(NT(dimension()));
        diam = 1.0 * diam;

        //std::cout<<"diam = "<<diam<<", radius = "<<radius<<std::endl;

    }

    template <class SpecSettings>
    void set_LMIatP_A0(SpecSettings& specSettings) {
        specSettings.LMIatP = getLMI().getA0();
    }

    void shift(VT e) {
        MT A0 = getLMI().getA0();
        std::vector<MT> matrices = getLMI().getMatrices();

        int d = matrices.size();

        for (int i = 0; i < d; ++i) {
            A0 = A0 + e(i)*matrices[i];
        }

        lmi.set_A0(A0);

    }

    void linear_transformIt(const MT &T){

        std::vector<MT> matrices;
        for (int i = 0; i < dimension(); ++i) {
            //lmi.evaluateWithoutA0_revised(T.col(i), LMIatP);
            lmi.evaluateWithoutA0(T.col(i), LMIatP);
            matrices.push_back(LMIatP);
        }
        for (int j = 0; j < dimension(); ++j) {
            lmi.set_Ai(matrices[j], j);
        }

    }

};

#endif //VOLESTI_SPECTRAHEDRON_H
