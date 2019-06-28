//
// Created by panagiotis on 26/6/2019.
//

#ifndef VOLESTI_SDP_PROBLEM_H
#define VOLESTI_SDP_PROBLEM_H

#include "Eigen"
#include "spectrahedron.h"
#include <vector>
#include <spectrahedron.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;

namespace optimization {

    class sdp_problem {
    private:
        Spectrahedron spectrahedron;
        VT objectiveFunction;
        VT solution;

    public:

        sdp_problem(std::istream &is) {

            std::string line;
            std::string::size_type sz;

            if (std::getline(is, line, '\n').eof())
                throw 1;

           if (!(!line.compare("minimize") || !line.compare("Minimize") || !line.compare("MINIMIZE") ||
                     !line.compare("MIN") || !line.compare("Min") || !line.compare("min")))
                throw 1;

            if (std::getline(is, line, '\n').eof())
                throw 1;

            std::list<double > obj;

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

            if (line.compare("subject to") && line.compare("SUBJECT TO") & line.compare("Subject To") &&
                line.compare("s.t."))
                throw 1;


            int matrixDim = -1;

            // for each matrix
            for (int atMatrix = 1; atMatrix <= matricesNum; atMatrix++) {
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


                    if (matrixDim == -1)
                        matrixDim = A.size();
                    else if (matrixDim != A.size())
                        throw  1;

                    rows.push_back(A);

                    //read first row of b
                } /*  while (!std::getline(is, line, '\n').eof())  */

                if (rows.size() != matrixDim)
                    throw  1;

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

            }

            LMI lmi(matrices);
            this->spectrahedron = Spectrahedron(lmi);
        }
    };
}
#endif //VOLESTI_SDP_PROBLEM_H
