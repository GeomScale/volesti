#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <chrono>
#include "Eigen/Eigen"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

// Assume that A is already symmetric

// Using LDLT decomposition: more numerically stable for singular matrices
bool isPosSemidefinite(MT A){
    Eigen::LDLT<MT> A_ldlt(A);
    if (A_ldlt.info() != Eigen::NumericalIssue && A_ldlt.isPositive())
        return true;
    return false;
}

// Using LLT decomposition
bool isPosSemidefinite2(MT A){
    Eigen::LLT<MT> A_llt(A);
    if(A_llt.info() != Eigen::NumericalIssue) return true;
    return false;
}

// Test 1: Reject Sampling

// Input:   n   : the matrix size
//          c   : a counter variable to keep the number of trials for the test  
MT AcceptReject(int n, int * c){
    int i,j;
    double val;
    MT A = MT::Identity(n,n);
    *c = 0;
    while(true){
        (*c)++;
        for(i = 0; i < n; i++){
            for(j = i+1; j < n; j++){
                val = 2* ((double) rand()/RAND_MAX) - 1;
                A(i,j) = val;
                A(j,i) = val;
            }
        }
        if(isPosSemidefinite(A)) return A;
    }
}

std::pair<int, double> testRejectSampling(int n, int num_matrices = 100){
    int numtrials = 0, c = 0;
    std::pair<int, double> ptest;
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < num_matrices; ++i){
        MT A = AcceptReject(n, &c);
        numtrials += c;
    }
    auto end = std::chrono::steady_clock::now();
    double timetest = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    ptest.first = numtrials;
    ptest.second = timetest;
    std::cout << "Matrix size : " << n << " : " << numtrials << " trials : " << timetest << " (ms)" << std::endl;
    std::cout << "===================" << std::endl;
    return ptest;
}

int main(){
    testRejectSampling(7);
    return 0;
}
The only difference in the case of correlation matrices is the per-step complexity of the random walk. 
Notice that there is no need for LMI evaluations, you could only store the current Markov point as a correlation matrix.
Similarly, with the direction HnR or Billiard Walk (or momenta in HMC), it would be a hollow matrix. 
Also with the computation of the gradient. Since matrices A_i are well-structured and super sparse,
the evaluation of the gradient is cheaper.

Please check this paper for details and try to specialize the per-step computation for the case of correlation matrices.
No need to change your proposal (it is pretty good for now), you could do this during the bonding period.



#include "matrix_operations/EigenvaluesProblems.h"


template <typename MT>
struct evaluate_lmi {
    
};

template <typename NT>
struct evaluate_lmi<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> > {
public:
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    MT vectorMatrix;

    int _m, _d;

    /// Create the vectorMatrix, which has at each column the distinct elements of each A_i, i=1,...,d
    void setVectorMatrix(int const& m, int const& d, std::vector<MT> &matrices) 
    {
        _m = m;
        _d = d;
        int newM = m * (m + 1) / 2;

        // allocate memory
        vectorMatrix.setZero(newM, d);
        //a.setZero(m);

        // initialze iterator and skip A_0
        typename std::vector<MT>::iterator iter = matrices.begin();
        iter++;

        // copy elements
        int atMatrix = 0;

        for (; iter != matrices.end(); iter++, atMatrix++) {
            int i = 0;

            for (int at_row = 0; at_row < m; at_row++)
                for (int at_col = at_row; at_col < m; at_col++) {
                    vectorMatrix(i++, atMatrix) = (*iter)(at_row, at_col);
                }

        }
    }