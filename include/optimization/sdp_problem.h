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
#include "lmi_strict_feasibility.h"
#include "interior_point_sdp.h"
#include "SDPA_format_manager.h"

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

        sdp_problem(Spectrahedron& spectrahedron, VT& objectiveFunction) {
            this->spectrahedron = spectrahedron;
            this->objectiveFunction = objectiveFunction;
        }

        sdp_problem(std::istream &is, bool SDPAformat = false, bool dual = false) {//TODO dual

            if (SDPAformat) {
                LMI lmi;
                SDPAFormat::loadSDPAFormatFile(is, lmi, objectiveFunction);
                this->spectrahedron = Spectrahedron(lmi);
                return;
            }

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

        void writeSDPAFormatFile(std::ostream& os) {
            LMI lmi;
            lmi = spectrahedron.getLMI();
            SDPAFormat::writeSDPAFormatFile(os, lmi, objectiveFunction);
        }


        VT getStrictlyFeasiblePoint() {
            VT point = getInteriorPoint(spectrahedron);//strict_feasible_point(spectrahedron, 0.00001, 0.00001, 100000);
            return point;
        }

        bool isStrictlyFeasible(VT& point) {
            LMI lmi;
            lmi = spectrahedron.getLMI();
            return lmi.isNegativeDefinite(point);
        }

        bool isFeasible(VT& point) {
            LMI lmi;
            lmi = spectrahedron.getLMI();
            return lmi.isNegativeSemidefinite(point);
        }

        VT getSolution() {
            return solution.first.getCoefficients();
        }

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
        void solve(Parameters &parameters, double error, unsigned int maxSteps, bool sampled_covariance_matrix=false) {
            Point initial(getStrictlyFeasiblePoint());
            if (!sampled_covariance_matrix)
                solution = cutting_plane_method(spectrahedron, objectiveFunction, parameters, error, maxSteps, initial);
            else
                solution = cutting_plane_method_sampled_covariance_matrix(spectrahedron, objectiveFunction, parameters, error, maxSteps, initial);
        }

        template <class Parameters>
        void solve(Parameters &parameters, double error, unsigned int maxSteps, Point& initial, bool sampled_covariance_matrix=false) {
            if (!sampled_covariance_matrix)
                solution = cutting_plane_method(spectrahedron, objectiveFunction, parameters, error, maxSteps, initial);
            else
                solution = cutting_plane_method_sampled_covariance_matrix(spectrahedron, objectiveFunction, parameters, error, maxSteps, initial);
        }

        void printSolution() {
            std::cout << "Min: " << solution.second << "\n";
            std::cout << "Coordinates:\n" << solution.first.getCoefficients() << "\n";
        }

        void saveToFile(std::ofstream& os) {
            os << "Minimize\n";

            for (int i=0 ; i<objectiveFunction.rows() ; i++)
                os << objectiveFunction(i) << " ";

            os << "\nSubject to\n";

            const LMI& lmi = spectrahedron.getLMI();
            const MT& A = lmi.getA0();

            for (int i=0 ; i<A.rows() ; i++) {
                for (int j=0 ; j<A.cols() ; j++)
                    os << A(i, j) << " ";

                os << "\n";
            }

            const std::vector<MT>& matrices = lmi.getMatrices();

            for (int matrix=0 ; matrix<matrices.size() ; matrix++) {
                MT A = matrices[matrix];

                for (int i=0 ; i<A.rows() ; i++) {
                    for (int j=0 ; j<A.cols() ; j++)
                        os << A(i, j) << " ";

                    os << "\n";
                }
            }
        }
    };
}
#endif //VOLESTI_SDP_PROBLEM_H
