
// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "samplers.h"
#include "rounding.h"
#include "sample_only.h"
#include "sdp_generator.h"
#include "spectrahedron.h"


typedef std::string::iterator string_it;
typedef std::list<double> listVector;


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


template <typename MT, typename LMII, typename VT>
void loadSDPAFormatFile(std::istream &is, LMII &lmi, VT &objectiveFunction) {
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

    while  (constantVector.size() < variablesNum) {
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

    lmi = LMII(matrices);
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector sdp_approx(Rcpp::Nullable<Rcpp::CharacterVector> file = R_NilValue,
               Rcpp::Nullable<int> N = R_NilValue,
                             Rcpp::Nullable<int> random_walk = R_NilValue,
                             Rcpp::Nullable<unsigned int> walk_length = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT, NT, VT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <MT, VT> lmi;
    typedef Spectrahedron <lmi, Point> spectaedro;

    std::ifstream inp;

    //std::string bar = "_";
    //std::string txt = ".txt";
    //std::string sdp = "sdp_prob"+bar+std::to_string(Rcpp::as<int>(d))+bar+std::to_string(Rcpp::as<int>(mm))+txt;

    //std::cout<<"reading... "<<sdp<<std::endl;

    inp.open(Rcpp::as<std::string>(file),std::ifstream::in);
    lmi Slmi;
    VT c;
    loadSDPAFormatFile<MT>(inp, Slmi, c);

    spectaedro SP(Slmi);//, SP2;
    unsigned int n = SP.dimension();
    //SP = generateSDP2<lmi, spectaedro, Point>(Rcpp::as<int>(nn), Rcpp::as<int>(mm));

    //SP2 = SP;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    typedef boost::mt19937 RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    bool round = false;
    std::pair <Point, NT> InnerB;
    Point p(n);
    NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0;
    InnerB.first = p;// = SP.ComputeInnerBall(diam_spec);

    vars <NT, RNGType> var(0, n, 1, 1, 0.0, 0.1, 0, 0.0, 0, InnerB.second, diam_spec, 0, false, rng, urdist,
                           urdist1,
                           -1.0, true, false, round, false, false, false, false, false, true);

    std::list <Point> randPoints, randPoints2;
    spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
    settings.LMIatP = SP.getLMI().getA0();
    preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
    settings.LMIatP = SP.getLMI().getA0();
    p = Point(n);

    spectaedro::BoundaryOracleBoltzmannHMCSettings settings2;
    settings2.first = true;
    settings2.epsilon = 0.0001;
    //settings2.LMIatP = SP.getLMI().getA0();
    Point cc(c);// = get_direction<RNGType, Point, NT>(n);

    //SP.print();
    //std::cout<<"n = ="<<n<<" c = "<<c<<std::endl;

    //std::filebuf fb;
    //fb.open("sdp_prob.txt", std::ios::out);
    //std::ostream os(&fb);
    //writeSDPAFormatFile<MT>(os, SP.getLMI(), c.get_coefficients());



    NT T = var.diameter;

    NT tred = 1.0 - 1.0/std::sqrt(NT(n));

    int walkL = 1, NN=50;
    if(walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);
    if(N.isNotNull()) NN = Rcpp::as<unsigned int>(N);
    //if(max_iterations.isNotNull()) max_iter = Rcpp::as<unsigned int>(max_iterations);
    Rcpp::NumericVector retvec(NN);

    if(!random_walk.isNotNull() || Rcpp::as<std::string>(random_walk).compare(std::string("HMC"))==0){


        //std::cout << "HMC"<< std::endl;
        //std::cout << optimal_solution<< std::endl;
        //s/td::cout << err<< std::endl;
        //std::cout << (std::abs(min_val - optimal_solution) / std::abs(optimal_solution)) << std::endl;

        for (int i = 0; i < NN; ++i) {
            for (int j = 0; j < walkL; ++j) {
                HMC_boltzmann_reflections(SP, p, diam_spec, var, cc, T, settings2);
            }
            retvec(i)= c.dot(p.getCoefficients());

            T = T * tred;

        }
    } else if(Rcpp::as<std::string>(random_walk).compare(std::string("RDHR"))==0) {


        //std::cout << "Hit and Run" << std::endl;

        for (int i = 0; i < NN; ++i) {
            for (int j = 0; j < walkL; ++j) {
                hit_and_run_Boltzmann_spec(p, SP, var, cc, T);
            }
            retvec(i)= c.dot(p.getCoefficients());
            T = T * tred;

        }
    }

    //std::cout<<"approx_value = "<<min_val<<std::endl;
    return retvec;

}