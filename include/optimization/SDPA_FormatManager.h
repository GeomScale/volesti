//
// Created by panagiotis on 20/7/2019.
//

#ifndef VOLESTI_SDPA_FORMAT_MANAGER_H
#define VOLESTI_SDPA_FORMAT_MANAGER_H

#include "spectrahedron.h"

#include <string>
#include <sstream>


typedef std::string::iterator string_it;
typedef std::list<double> listVector;


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


    template<typename NT, typename MT, typename VT>
    void loadSDPAFormatFile(std::ifstream &is, LMI<NT, MT, VT> &lmi, VT &objectiveFunction) {
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

        while (constantVector.size() < variablesNum) {
            if (std::getline(is, line, '\n').eof())
                throw 1;
            listVector t = readVector(line);
            constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
        }

//        for (auto x : constantVector)
//            std::cout << x << "  ";
//        std::cout << "\n";
//            throw 1;


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
//                            std::cout <<value << " val\n";
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

        lmi = LMI<NT, MT, VT>(matrices);
    }




#endif //VOLESTI_SDPA_FORMAT_MANAGER_H
