//
// Created by panagiotis on 27/5/2019.
//

#ifndef VOLESTI_LP_PROBLEM_H
#define VOLESTI_LP_PROBLEM_H

#include "polytopes.h"
#include "Eigen"
#include "interior_point.h"

namespace optimization {

    typedef enum Goal {
        minimize, maximize
    } Goal;

    template<class Point, typename NT>
    class lp_problem {
    public:
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

        HPolytope<Point> polytope;
        VT objectiveFunction;
        Goal goal;
        Point solution;
        NT solutionVal;

        lp_problem(HPolytope<Point>& hPolytope, VT& objectiveFunction, Goal goal) {
            this->polytope = hPolytope;
            this->objectiveFunction = objectiveFunction;
            this->goal = goal; //TODO now it only minimuzes
        }

        lp_problem() {}

        // CPLEX LP file format
        lp_problem(std::istream &is) {

            std::string line;
            std::string::size_type sz;

            if (std::getline(is, line, '\n').eof())
                throw 1;

            if (!line.compare("maximize") || !line.compare("Maximize") || !line.compare("MAXIMIZE") || !line.compare("MAX") || !line.compare("Max") || !line.compare("max"))
                goal = maximize;
            else if (!line.compare("minimize") || !line.compare("Minimize") || !line.compare("MINIMIZE") || !line.compare("MIN") || !line.compare("Min") || !line.compare("min"))
                goal = minimize;
            else
                throw 1;

            if (std::getline(is, line, '\n').eof())
                throw 1;

            std::list<NT> obj;

            try {
                NT num = std::stod(line, &sz);
                obj.push_back(num);

                while (1) {
                    line = line.substr(sz);
                    num = std::stod(line, &sz);
                    obj.push_back(num);
                }
            }
            catch (std::exception& e) {
                if (obj.empty())
                    throw 1;
            }

            unsigned int dim = obj.size();
            this->objectiveFunction.resize(dim);

            int i=0;
            for (auto x : obj) {
                this->objectiveFunction(i) = x;
                i++;
            }

            if (std::getline(is, line, '\n').eof())
                throw 1;

            if (line.compare("subject to") && line.compare("SUBJECT TO") & line.compare("Subject To") && line.compare("s.t."))
                throw 1;

            //read constraints
            std::list<std::list<NT> > constraints;
            std::list<NT> b;

            // for each constraint
            while (!std::getline(is, line, '\n').eof()) {
                std::list<NT> A;
                // read first line of A

                NT num = std::stod(line, &sz);
                A.push_back(num);

                for (int j=2 ; j<=dim  ; j++) {
                    line = line.substr(sz);
                    num = std::stod(line, &sz);
                    A.push_back(num);
                }

                constraints.push_back(A);

                //read first row of b
                line = line.substr(sz);
                num = std::stod(line, &sz);
                b.push_back(num);
            }

            MT A;
            VT _b;
            A.resize(constraints.size(), dim);
            _b.resize(constraints.size());

            i=0;
            for (auto constraint : constraints) {
                int j = 0;

                for (auto coeff : constraint) {
                    A(i, j) = coeff;
                    j++;
                }
                i++;
            }

            i=0;
            for (auto coeff : b) {
                _b(i) = coeff;
                i++;
            }

            this->polytope.init(dim, A, _b);
        }


        template<class Parameters>
        void solve(Parameters parameters, const NT error, unsigned int maxSteps, bool sampledCovarianceMatrix = false) {
            Point initial = getInteriorPoint<Point, NT>(polytope);
            std::pair<Point, NT> sol;


            if (sampledCovarianceMatrix)
                sol = cutting_plane_method_sampled_covariance_matrix(polytope, objectiveFunction, parameters, error, maxSteps,
                                                          initial);
            else
                sol = cutting_plane_method(polytope, objectiveFunction, parameters, error, maxSteps, initial);

            solutionVal = sol.second;
            solution = sol.first;
        }


        const HPolytope<Point> &getHPolytope() const {
            return polytope;
        }

        void setHPolytope(const HPolytope<Point> &hPolytope) {
            lp_problem::polytope = hPolytope;
        }

        const VT &getObjectiveFunction() const {
            return objectiveFunction;
        }

        void setObjectiveFunction(const VT &objectiveFunction) {
            lp_problem::objectiveFunction = objectiveFunction;
        }

        Goal getGoal() const {
            return goal;
        }

        void setGoal(Goal goal) {
            lp_problem::goal = goal;
        }

        const HPolytope<Point> &getPolytope() const {
            return polytope;
        }

        void setPolytope(const HPolytope<Point> &polytope) {
            lp_problem::polytope = polytope;
        }

        Point getSolution() const {
            return solution;
        }

        void setSolution(Point solution) {
            lp_problem::solution = solution;
        }

        NT getSolutionVal() const {
            return solutionVal;
        }

        void setSolutionVal(NT solutionVal) {
            lp_problem::solutionVal = solutionVal;
        }

        void printSolution() {
            std::cout << "Min value is: " << solutionVal << std::endl <<
                      "coords: ";
            solution.print();
        }

        void saveToFile(std::ofstream& os) {
            if (goal == minimize)
                os << "Minimize\n";
            else
                os << "Maximize\n";

            for (int i=0 ; i<objectiveFunction.rows() ; i++)
                os << objectiveFunction(i) << " ";

            os << "\nSubject to\n";

            MT A = polytope.get_mat();
            VT b = polytope.get_vec();

            for (int i=0 ; i<A.rows() ; i++) {
                for (int j=0 ; j<A.cols() ; j++)
                    os << A(i, j) << " ";

                os << b(i) << "\n";
            }

        }
    };
}

#endif //VOLESTI_LP_PROBLEM_H
