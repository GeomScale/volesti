#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <chrono>
#include "Eigen/Eigen"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;

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
