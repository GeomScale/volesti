#ifndef SPEED_TESTS_H
#define SPEED_TESTS_H

#include "ChebTools/ChebTools.h"
#include <map>

double real_roots_time(ChebTools::ChebyshevExpansion &ce, long N);
double plus_by_inplace(ChebTools::ChebyshevExpansion &ce, const ChebTools::ChebyshevExpansion &ce2, int N);
double mult_by_inplace(ChebTools::ChebyshevExpansion &ce, double val, int N);
void mult_by(ChebTools::ChebyshevExpansion &ce, double val, int N);
std::map<std::string, double> evaluation_speed_test(ChebTools::ChebyshevExpansion &cee, const Eigen::VectorXd &xpts, long N) ;
Eigen::MatrixXd eigs_speed_test(std::vector<std::size_t> &Nvec, std::size_t Nrepeats);

#endif