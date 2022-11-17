//
// Created by panagiotis on 9/7/2019.
//

#ifndef VOLESTI_SDP_GENERATOR_H
#define VOLESTI_SDP_GENERATOR_H

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

typedef boost::mt19937 RNGType;

/// Generates a random matrix
/// @tparam NT Number type
template <class NT>
void randomMatrixGOE(Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>& M) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    unsigned m = M.rows();
    boost::normal_distribution<> rdist(0,1);
    unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();//4 if fixed for debug
    RNGType rng(seed);

    for (unsigned int i=0; i<m; i++) {
        for (unsigned int j=0; j<m; j++) {
            M(i,j) = rdist(rng);
        }
    }
}

/// Generates a random spectahedron S(n, m)
/// @tparam NT Number type
template<typename NT>
Spectrahedron<NT, Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<NT, Eigen::Dynamic, 1> >  generateSDP(int n, int m) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    MT ones = MT::Ones(m, m);
    MT M = 2* Eigen::MatrixXd::Random(m,m) - ones;

    MT I = Eigen::MatrixXd::Identity(m, m);
    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;

    std::cout<<"A0 = "<<matrices[0]<<"\n"<<std::endl;


    MT ones2 = MT::Ones(m/2, m/2);
    MT MM(m/2, m/2), MMM(m/2, m/2);

    for (int i=1 ; i<=n ; i++) {
        MM =  2 * MT::Random(m/2, m/2) - ones2;
        std::cout<<"MM = "<<MM<<"\n"<<std::endl;
        std::cout<<"MM.transpose() = "<<MM.transpose()<<"\n"<<std::endl;

        MMM.setZero(m/2, m/2);
        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m / 2; ++k) {
                MMM(j,k) = MM(j,k) + MM(k,j);
            }
        }

        //MM = MM + MM.transpose();

        std::cout<<"MMM = "<<MMM<<"\n"<<std::endl;

        MT A;
        A.setZero(m, m);

        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m/2; ++k) {
                A(j,k) = MMM(j,k);
                A(j+m/2, k+m/2) = -MMM(j, k);
            }
        }
        std::cout<<"A = "<<A<<"\n"<<std::endl;
        matrices[i] = A;
    }

    LMI<NT, MT, VT> lmi(matrices);
    Spectrahedron<NT, MT, VT> spectrahedron(lmi);
    return spectrahedron;

    //return optimization::sdp_problem<Point>(spectrahedron, obj);
}

/// Generates a random spectahedron S(n, m)
/// @tparam NT Number type
template<typename NT>
Spectrahedron<NT, Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<NT, Eigen::Dynamic, 1> >  generateSDP2(int n, int m) {

    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    MT ones = MT::Ones(m, m);
    MT M = 2* Eigen::MatrixXd::Random(m,m) - ones;

    MT I = Eigen::MatrixXd::Identity(m, m);
    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;

    //std::cout<<"A0 = "<<matrices[0]<<"\n"<<std::endl;


    MT ones2 = MT::Ones(m/2, m/2);
    MT MM(m/2, m/2), MMM(m/2, m/2);

    for (int i=1 ; i<=n ; i++) {
        MM =  2 * MT::Random(m/2, m/2) - ones2;
        MMM.setZero(m/2, m/2);
        for (int l = 0; l < m/2; ++l) {
            MMM(l,l) = MM(l,l);
        }
        //std::cout<<"MM = "<<MM<<"\n"<<std::endl;
        //std::cout<<"MM.transpose() = "<<MM.transpose()<<"\n"<<std::endl;


        for (int j = 0; j < m/2; ++j) {
            for (int k = j+1; k < m / 2; ++k) {
                MMM(j,k) = MM(j,k);
                MMM(k,j) = MMM(j,k);
            }
        }

        //MM = MM + MM.transpose();

        //std::cout<<"MMM = "<<MMM<<"\n"<<std::endl;

        MT A;
        A.setZero(m, m);

        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m/2; ++k) {
                A(j,k) = MMM(j,k);
                A(j+m/2, k+m/2) = -MMM(j, k);
            }
        }
        //std::cout<<"A = "<<A<<"\n"<<std::endl;
        matrices[i] = A;
    }

    LMI<NT, MT, VT> lmi(matrices);
    Spectrahedron<NT, MT, VT> spectrahedron(lmi);
    return spectrahedron;

    //return optimization::sdp_problem<Point>(spectrahedron, obj);
}


#endif //VOLESTI_SDP_GENERATOR_H
