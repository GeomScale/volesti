// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SDPA_FORMAT_MANAGER_H
#define VOLESTI_SDPA_FORMAT_MANAGER_H

#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "convex_bodies/spectrahedra/sparse_spectrahedron.h"

#include <string>
#include <sstream>
#include <iostream>

/// Reads/writes files according to the SDPA format for SDPs.
/// Supports both DENSE and SPARSE formats:
///
/// DENSE FORMAT:
/// <number of variables>
/// <number of blocks>
/// <block structure>
/// <objective function vector>
/// <all entries of each matrix, row by row>
///
/// SPARSE FORMAT (SDPA standard):
/// <number of variables>
/// <number of blocks>
/// <block structure>
/// <objective function vector>
/// <matno> <blkno> <i> <j> <value>  (only non-zero entries)
///
/// \tparam NT Numerical Type
template <typename NT>
class SdpaFormatManager {
private:
    typedef std::string::iterator string_it;
    typedef std::list<NT> listVector;

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    /// Return the first non white space/tab character and advance the iterator one position
    /// @param[in, out] it Current position
    /// @param[in] end End of string
    /// @return First non white space/tab character
    char consumeSymbol(string_it &at, string_it & end) {
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

    /// Determine if current line is a comment
    /// @param[in] line The current line
    /// @return true if line is a comment, false otherwise
    bool isCommentLine(std::string & line) {
        string_it at = line.begin();
        string_it end = line.end();
        char c = consumeSymbol(at, end);
        return c == '"' || c == '*' || c == '#' || line.empty();
    }

    /// Get an integer from the string
    /// \param[in] string
    /// \return an integer
    int fetchNumber(std::string &string) {
        std::stringstream stream(string);
        int num;
        stream >> num;
        return num;
    }

    /// Read a vector of the form {val1, val2, ..., valn}
    /// @param string Contains the vector
    /// @return a list with the n numbers
    listVector readVector(std::string &string) {
        std::stringstream stream(string);
        listVector vector;
        NT value;

        while (stream >> value) {
            vector.push_back(value);
        }

        return vector;
    }

    /// Detect if file is in sparse or dense format
    /// Sparse format has entries like: "matno blkno i j value"
    /// Dense format has just numerical values
    /// @param[in] is Input stream
    /// @param[out] firstDataLine First data line after header
    /// @return true if sparse format, false if dense
    bool detectSparseFormat(std::ifstream &is, std::string &firstDataLine) {
        // Save current position
        std::streampos start_pos = is.tellg();
        
        std::string line;
        
        // Skip header lines (comments, variables, blocks, structure, objective)
        int lines_to_skip = 0;
        while (std::getline(is, line)) {
            if (!isCommentLine(line)) {
                lines_to_skip++;
                if (lines_to_skip >= 4) break; // After objective function
            }
        }
        
        // Read first data line
        while (std::getline(is, line)) {
            if (!line.empty() && !isCommentLine(line)) {
                firstDataLine = line;
                break;
            }
        }
        
        // Restore file position
        is.clear();
        is.seekg(start_pos);
        
        // Check format of first data line
        std::istringstream iss(firstDataLine);
        int count = 0;
        NT value;
        while (iss >> value) {
            count++;
        }
        
        // Sparse format: exactly 5 values (matno blkno i j value)
        // Dense format: many values or just a few (matrix entries)
        // Heuristic: if exactly 5 integers, likely sparse
        iss.clear();
        iss.str(firstDataLine);
        int int_count = 0;
        int temp;
        while (iss >> temp && int_count < 4) {
            int_count++;
        }
        
        bool is_sparse = (count == 5 && int_count == 4);
        
        if (is_sparse) {
            std::cout << "SPARSE SDPA format" << std::endl;
        } else {
            std::cout << "DENSE SDPA format" << std::endl;
        }
        
        return is_sparse;
    }

public:

    /// Reads an SDPA format file in DENSE format
    /// \param[in] is An open stream pointing to the file
    /// \param[out] matrices the matrices A0, A1, A2, ..., An
    /// \param[out] objectiveFunction The objective function of the sdp
    void loadSDPAFormatFileDense(std::ifstream &is, std::vector<MT> &matrices, VT &objectiveFunction) {
        std::string line;

        std::getline(is, line, '\n');

        //skip comments
        while (isCommentLine(line)) {
            std::getline(is, line, '\n');
        }

        //read variables number
        int variablesNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read number of blocks
        int blocksNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read block structure vector
        listVector blockStructure = readVector(line);

        if (blockStructure.size() != blocksNum)
            throw std::runtime_error("Wrong number of blocks");

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read constant vector (objective function)
        listVector constantVector = readVector(line);

        while (constantVector.size() < variablesNum) {
            if (std::getline(is, line, '\n').eof())
                throw std::runtime_error("Unexpected end of file");

            listVector t = readVector(line);
            constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
        }

        matrices = std::vector<MT>(variablesNum + 1);
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
                            throw std::runtime_error("Unexpected end of file while reading matrix");

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
                    }

                } else { //read diagonal block
                    blockSize = std::abs(blockSize);
                    int at = 0;

                    while (at < blockSize) {
                        if (std::getline(is, line, '\n').eof())
                            throw std::runtime_error("Unexpected end of file while reading diagonal block");

                        listVector vec = readVector(line);

                        for (double value : vec) {
                            matrix(offset + at, offset + at) = value;
                            at++;
                        }
                    }
                }

                offset += std::abs(blockSize);
            }

            //the LMI in SDPA format is >0, I want it <0
            if (atMatrix == 0) //F0 has - before it in SDPA format, the rest have +
                matrices[atMatrix] = matrix;
            else
                matrices[atMatrix] = -1 * matrix;
        }

        // return objective function
        objectiveFunction.setZero(variablesNum);
        int at = 0;

        for (auto value : constantVector)
            objectiveFunction(at++) = value;
    }

    /// Reads an SDPA format file in SPARSE format
    /// Format: <matno> <blkno> <i> <j> <value>
    /// \param[in] is An open stream pointing to the file
    /// \param[out] matrices the matrices A0, A1, A2, ..., An
    /// \param[out] objectiveFunction The objective function of the sdp
    void loadSDPAFormatFileSparse(std::ifstream &is, std::vector<MT> &matrices, VT &objectiveFunction) {
        std::string line;

        std::getline(is, line, '\n');

        //skip comments
        while (isCommentLine(line)) {
            std::getline(is, line, '\n');
        }

        //read variables number
        int variablesNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read number of blocks
        int blocksNum = fetchNumber(line);

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read block structure vector
        listVector blockStructure = readVector(line);

        if (blockStructure.size() != blocksNum)
            throw std::runtime_error("Wrong number of blocks");

        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        //read objective function
        listVector constantVector = readVector(line);

        while (constantVector.size() < variablesNum) {
            if (std::getline(is, line, '\n').eof())
                throw std::runtime_error("Unexpected end of file");

            listVector t = readVector(line);
            constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
        }

        // Calculate matrix dimension
        int matrixDim = 0;
        for (auto x : blockStructure)
            matrixDim += std::abs((int) x);

        // Initialize matrices with zeros
        matrices = std::vector<MT>(variablesNum + 1);
        for (int i = 0; i <= variablesNum; ++i) {
            matrices[i].setZero(matrixDim, matrixDim);
        }

        // Read sparse entries
        // Format: matno blkno i j value
        // matno: 0 for A0, 1..n for A1..An
        // blkno: block number (1-indexed, but we only support 1 block)
        // i, j: row, column (1-indexed)
        // value: matrix entry
        
        int entries_read = 0;
        while (std::getline(is, line)) {
            if (line.empty() || isCommentLine(line)) continue;
            
            std::istringstream iss(line);
            int matno, blkno, i, j;
            NT value;
            
            if (iss >> matno >> blkno >> i >> j >> value) {
                // Convert from 1-indexed to 0-indexed
                i--; 
                j--;
                
                // Validate indices
                if (matno < 0 || matno > variablesNum) {
                    std::cerr << "Warning: Invalid matrix number " << matno << std::endl;
                    continue;
                }
                if (i < 0 || i >= matrixDim || j < 0 || j >= matrixDim) {
                    std::cerr << "Warning: Invalid indices (" << i << "," << j 
                              << ") for matrix size " << matrixDim << std::endl;
                    continue;
                }
                
                // Store entry with correct sign convention
                // SDPA format: constraint is F0 - sum(xi * Fi) >= 0
                // Our format: A0 + sum(xi * Ai) <= 0
                // So: A0 = -F0, Ai = Fi
                if (matno == 0) {
                    matrices[0](i, j) = value;
                    if (i != j) {
                        matrices[0](j, i) = value; // Symmetric
                    }
                } else {
                    matrices[matno](i, j) = -value;
                    if (i != j) {
                        matrices[matno](j, i) = -value; // Symmetric
                    }
                }
                
                entries_read++;
            }
        }
        
        std::cout << "Read " << entries_read << " sparse entries" << std::endl;

        // Set objective function
        objectiveFunction.setZero(variablesNum);
        int at = 0;
        for (auto value : constantVector)
            objectiveFunction(at++) = value;
    }

    /// Reads an SDPA format file (auto-detects dense or sparse)
    /// \param[in] is An open stream pointing to the file
    /// \param[out] matrices the matrices A0, A1, A2, ..., An
    /// \param[out] objectiveFunction The objective function of the sdp
    void loadSDPAFormatFile(std::ifstream &is, std::vector<MT> &matrices, VT &objectiveFunction) {
        // Save position to restart
        std::streampos start = is.tellg();
        
        // Try to detect format
        std::string firstDataLine;
        bool is_sparse = detectSparseFormat(is, firstDataLine);
        
        // Reset to beginning
        is.clear();
        is.seekg(start);
        
        // Load with appropriate method
        if (is_sparse) {
            loadSDPAFormatFileSparse(is, matrices, objectiveFunction);
        } else {
            loadSDPAFormatFileDense(is, matrices, objectiveFunction);
        }
    }

    /// Create a SDPA format file (dense format)
    /// \param[in] os Open stream to file
    /// \param[in] matrices The matrices A0, ..., An
    /// \param[in] objectiveFunction The objective function of the sdp
    void writeSDPAFormatFile(std::ostream &os, std::vector<MT> const & matrices, VT const & objectiveFunction) {
        int dim = matrices.size() - 1;
        MT A0 = matrices[0];

        os << dim << "\n";
        os << 1 << "\n";
        os << A0.rows() << "\n";

        os << objectiveFunction.transpose() << "\n";

        for (int i = 0; i < A0.rows(); i++)
            os << A0.row(i) << "\n";

        for (int at=1 ; at<matrices.size() ; ++at)
            for (int i = 0; i < matrices[at].rows(); i++)
                os << -1 * matrices[at].row(i) << "\n";
    }

    /// Read a spectrahedron and a vector (objective function) from a SDPA format input file
    /// \tparam Point
    /// \param[in] is opened stream to input file
    /// \param[out] spectrahedron
    /// \param[out] objectiveFunction
    template <typename Spectrahedron, typename Point>
    void loadSDPAFormatFile(std::ifstream &is, Spectrahedron &spectrahedron, Point &objectiveFunction) {
        std::vector<MT> matrices;
        VT coeffs;
        loadSDPAFormatFile(is, matrices, coeffs);
        LMI<NT, MT, VT> lmi(matrices);
        spectrahedron = Spectrahedron(lmi);
        objectiveFunction = Point(coeffs);
    }

    /// Write a spectrahedron and a vector (objective function) to a SDPA format output file
    /// \tparam Point
    /// \param[in] os opened stream to output file
    /// \param[in] spectrahedron
    /// \param[in] objectiveFunction
    template <typename Spectrahedron, typename Point>
    void writeSDPAFormatFile(std::ostream &os, Spectrahedron const & spectrahedron, Point const & objectiveFunction) {
        writeSDPAFormatFile(os, spectrahedron.getLMI().getMatrices(), objectiveFunction.getCoefficients());
    }
};

#endif //VOLESTI_SDPA_FORMAT_MANAGER_H