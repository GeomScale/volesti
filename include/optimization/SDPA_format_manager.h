//
// Created by panagiotis on 20/7/2019.
//

#ifndef VOLESTI_SDPA_FORMAT_MANAGER_H
#define VOLESTI_SDPA_FORMAT_MANAGER_H

#include "spectrahedron.h"

#include "Eigen"
#include <string>
#include <sstream>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef std::string::iterator string_it;
typedef std::list<double> listVector;

namespace SDPAFormat {

/**
 * Return the first non white space/tab character and advance the iterator one position
 * @param it
 * @return
 */
    char consumeSymbol(string_it &at, string_it &end) {
        while (at != end) {
            if (*at != ' ' && *at != '\t') {
                char c = *at;
                at++;
                return c;
            }

            at++;
        }

        return '\0';
    }


    bool isCommentLine(std::string &line) {
        string_it at = line.begin();
        string_it end = line.end();

        char c = consumeSymbol(at, end);

        return c == '"' || c == '*';
    }


    int fetchNumber(std::string &string) {
        std::stringstream stream(string);
        int num;
        stream >> num;
        return num;
    }


/**
 * Read a vector of the form {val1, val2, ..., valn}
 * @param string
 * @return
 */
    listVector readVector(std::string &string) {
        std::stringstream stream(string);
        listVector vector;
        double value;

        while (stream >> value) {
            vector.push_back(value);
        }

        return vector;
    }


    void loadSDPAFormatFile(std::istream &is, LMI &lmi, VT &objectiveFunction) {
        std::string line;
        std::string::size_type sz;

        std::getline(is, line, '\n');

        //skip comments
        while (isCommentLine(line)) {
            std::getline(is, line, '\n');
        }

        //read variables number
        int variablesNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw 1;

        //read number of blocks
        int blocksNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw 1;

        //read block structure vector
        listVector blockStructure = readVector(line); //TODO different if we have one block

        if (blockStructure.size() != blocksNum)
            throw 1;

        if (std::getline(is, line, '\n').eof())
            throw 1;

        //read constant vector
        listVector constantVector = readVector(line);

        if (constantVector.size() != variablesNum)
            throw 1;


        std::vector<MT> matrices(variablesNum + 1);
        int matrixDim = 0;
        for (auto x : blockStructure)
            matrixDim += std::abs((int) x);

        //read constraint matrices
        for (int atMatrix = 0; atMatrix < matrices.size(); atMatrix++) {
            MT matrix;
            matrix.setZero(matrixDim, matrixDim);

            int offset = 0;

            for (auto blockSize : blockStructure) {

                if (blockSize > 0) { //read a block blockSize x blockSize
                    int at = 0;
                    int i = 0, j = 0;

                    while (at < blockSize * blockSize) {
                        if (std::getline(is, line, '\n').eof())
                            throw 1;

                        listVector vec = readVector(line);

                        for (double value : vec) {
                            matrix(offset + i, offset + j) = value;
                            at++;
                            if (at % (int) blockSize == 0) { // new row
                                i++;
                                j = 0;
                            } else { //new column
                                j++;
                            }
                        }
                    } /* while (at<blockSize*blockSize) */

                } else { //read diagonal block
                    blockSize = std::abs(blockSize);
                    int at = 0;

                    while (at < blockSize) {
                        if (std::getline(is, line, '\n').eof())
                            throw 1;

                        listVector vec = readVector(line);

                        for (double value : vec) {
                            matrix(offset + at, offset + at) = value;
                            at++;
                        }
                    } /* while (at<blockSize) */
                }

                offset += std::abs(blockSize);
            } /* for (auto blockSize : blockStructure) */

            //the LMI in SDPA format is >0, I want it <0
            if (atMatrix == 0) //F0 has - before it in SDPA format, the rest have +
                matrices[atMatrix] = matrix;
            else
                matrices[atMatrix] = -1 * matrix;
        }

        // return lmi and objective function
        objectiveFunction.setZero(variablesNum);
        int at = 0;

        for (auto value : constantVector)
            objectiveFunction(at++) = value;

        lmi = LMI(matrices);
    }

    void writeSDPAFormatFile(std::ostream &os, LMI &lmi, VT &objectiveFunction) {
        int dim = lmi.getDim();
        MT A0 = lmi.getA0();
        std::vector<MT> matrices = lmi.getMatrices();
        os << dim << "\n";
        os << 1 << "\n";
        os << A0.rows() << "\n";

        os << objectiveFunction.transpose() << "\n";

        for (int i = 0; i < A0.rows(); i++)
            os << A0.row(i) << "\n";

        for (MT matrix : matrices)
            for (int i = 0; i < matrix.rows(); i++)
                os << -1 * matrix.row(i) << "\n";
    }

}
#endif //VOLESTI_SDPA_FORMAT_MANAGER_H
