#include <iostream>
#include <Eigen/Dense>
#include "preprocess/chebyshev_highs.hpp"

int main() {

    Eigen::MatrixXd A(6,3);
    A << 1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1;

    Eigen::VectorXd b = Eigen::VectorXd::Ones(6);

    Eigen::VectorXd center;
    double r;

    bool ok = volesti::compute_chebyshev_highs(A,b,center,r);

    if(!ok) {
        std::cout << "Failed\n";
        return 1;
    }

    std::cout << "Center: " << center.transpose() << "\n";
    std::cout << "Radius: " << r << "\n";

    return 0;
}

// same cube example, call multiple solvers
// (HiGHS implementation + placeholder CLP if linked)