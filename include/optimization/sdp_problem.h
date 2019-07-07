//
// Created by panagiotis on 26/6/2019.
//

#ifndef VOLESTI_SDP_PROBLEM_H
#define VOLESTI_SDP_PROBLEM_H

#include "Eigen"
#include "spectrahedron.h"
#include <vector>
#include <spectrahedron.h>
#include "lp_problem.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;

namespace optimization {

    template <class Point>
    class sdp_problem {
    private:
        Spectrahedron spectrahedron;
        VT objectiveFunction;
        std::pair<Point, double> solution;

        void read_inequality_formulation(std::istream &is) {
            std::string line;
            std::string::size_type sz;

            if (std::getline(is, line, '\n').eof())
                throw 1;

            if (!(!line.compare("minimize") || !line.compare("Minimize") || !line.compare("MINIMIZE") ||
                  !line.compare("MIN") || !line.compare("Min") || !line.compare("min")))
                throw 1;

            if (std::getline(is, line, '\n').eof())
                throw 1;

            std::list<double> obj;

            try {
                double num = std::stod(line, &sz);
                obj.push_back(num);

                while (1) {
                    line = line.substr(sz);
                    num = std::stod(line, &sz);
                    obj.push_back(num);
                }
            }
            catch (std::exception &e) {
                if (obj.empty())
                    throw 1;
            }

            unsigned int dim = obj.size();
            this->objectiveFunction.resize(dim);

            int i = 0;
            for (auto x : obj) {
                this->objectiveFunction(i) = x;
                i++;
            }


            int matricesNum = dim + 1;
            std::vector<MT> matrices(matricesNum);

            ///////////////////////////////////////////////////////////////////////////////////
            ///// Now read #matricesNum square matrices

            if (std::getline(is, line, '\n').eof())
                throw 1;

            if (!(line.compare("subject to") && line.compare("SUBJECT TO") && line.compare("Subject To") &&
                  line.compare("s.t.")))
                throw 1;


            int matrixDim = -1;

            // for each matrix
            for (int atMatrix = 0; atMatrix < matricesNum; atMatrix++) {
                int atLine = 1;
                std::list<std::list<double> > rows;

                // read matrix
                while (!std::getline(is, line, '\n').eof()) {
                    std::list<double> A;

                    double num = std::stod(line, &sz);
                    A.push_back(num);

                    try {
                        while (1) {
                            line = line.substr(sz);
                            num = std::stod(line, &sz);
                            A.push_back(num);
                        }
                    }
                    catch (std::exception &e) {
                        if (A.empty())
                            throw 1;
                    }


                    // enough elements in row
                    if (matrixDim == -1)
                        matrixDim = A.size();
                    else if (matrixDim != A.size())
                        throw 1;

                    rows.push_back(A);

                    //completed matrix
                    if (rows.size() == matrixDim)
                        break;

                } /*  while (!std::getline(is, line, '\n').eof())  */

                if (rows.size() != matrixDim)
                    throw 1;

                MT A(matrixDim, matrixDim);

                i = 0;
                for (auto row : rows) {
                    int j = 0;

                    for (auto coeff : row) {
                        A(i, j) = coeff;
                        j++;
                    }
                    i++;
                }

                matrices[atMatrix] = A;
            }

            LMI lmi(matrices);
            this->spectrahedron = Spectrahedron(lmi);
        }

    public:

        sdp_problem() {}

        sdp_problem(std::istream &is, bool dual = false) {

            if (!dual) {
                read_inequality_formulation(is);
            } else {
                //TODO
            }
        }

        void print() {
            std::cout << "min\n" << this->objectiveFunction.transpose() << "\n s.t.\n";
            this->spectrahedron.print();
        }


//        template<class Parameters>
//        void feasible_pont(Point &point, Parameters &parameters, double error, unsigned int maxSteps) {
//            LMI lmi;
//            lmi = spectrahedron.getLMI();
//            MT A0 = lmi.getA0();
//            int m = A0.rows();
//
//            Eigen::EigenSolver<MT> solver;
//            solver.compute(A0);
//            Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
//
//            double max_gamma = eivals(0).real();//TODO check if complex / real
//            std::cout << eivals << "\n";
//            for (int i = 1; i < eivals.rows(); i++)
//                if (eivals(i).real() > max_gamma)
//                    max_gamma = eivals(i).real();
//
//            MT matrix;
//            matrix.setZero(m,m);
//            for (int i=0 ; i<m ; i++)
//                matrix(i,i) = -1;
//
//            lmi.addMatrix(matrix);
//            Spectrahedron s(lmi);
//            VT objectiveFunction;
//            objectiveFunction.setZero(this->objectiveFunction.rows() + 1);
//            objectiveFunction(this->objectiveFunction.rows()) = 1;
//
//            VT zero;
//            zero.setZero(this->objectiveFunction.rows() + 1);
//            zero(this->objectiveFunction.rows()) = max_gamma;
//            Point initial(zero);
//            std::cout << max_gamma << "\n";
//            std::pair<Point, double> ret = cutting_plane_method(s, objectiveFunction, parameters, error, maxSteps, initial);
//            std::cout << ret.first.getCoefficients() << "\n";
//            assert(ret.first.getCoefficients()(ret.first.getCoefficients().rows() -1) <= 0);
//
//            VT _point(ret.first.getCoefficients().rows() - 1);
//            for (int i=0 ; i<ret.first.getCoefficients().rows() - 1 ; i++)
//                _point(i) = ret.first.getCoefficients()(i);
//
//            point = Point(_point);
//        }

        void transformFromLP(std::istream &is) {
            lp_problem<Point, double > lp(is);
            HPolytope<Point> polytope = lp.getHPolytope();
            MT A = polytope.get_mat();
            VT b = polytope.get_vec();

            objectiveFunction = lp.objectiveFunction;

            int matricesNum = A.cols() + 1;
            int dim = b.rows();

            std::vector<MT> matrices(matricesNum);

            // create A0
            MT matrix;
            matrix.setZero(dim, dim);

            for (int j=0 ; j<dim ; j++)
                matrix(j,j) = -b(j);

            matrices[0] = matrix;

            for (int i=1 ; i<matricesNum ; i++) {
                matrix.setZero(dim, dim);
                VT a = A.col(i-1);

                for (int j=0 ; j<dim ; j++)
                    matrix(j,j) = a(j);

                matrices[i] = matrix;
            }

            LMI lmi(matrices);
            this->spectrahedron = Spectrahedron(lmi);
        }

        template <class Parameters>
        void solve(Parameters &parameters, double error, unsigned int maxSteps) {
            Point initial(objectiveFunction.rows());
            solution = cutting_plane_method(spectrahedron, objectiveFunction, parameters, error, maxSteps, initial);
        }

        void printSolution() {
            std::cout << "Min: " << solution.second << "\n";
            std::cout << "Coordinates: " << solution.first.getCoefficients().transpose() << "\n";
        }
    };
}
#endif //VOLESTI_SDP_PROBLEM_H
