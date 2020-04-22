// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
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


typedef std::string::iterator string_it3;
typedef std::list<double> listVector3;


char consumeSymbol3(string_it3 &at, string_it3 &end) {
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


bool isCommentLine3(std::string &line) {
    string_it3 at = line.begin();
    string_it3 end = line.end();

    char c = consumeSymbol3(at, end);

    return c == '"' || c == '*';
}


int fetchNumber3(std::string &string) {
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
listVector3 readVector3(std::string &string) {
    std::stringstream stream(string);
    listVector3 vector;
    double value;

    while (stream >> value) {
        vector.push_back(value);
    }

    return vector;
}


template <typename MT, typename LMII, typename VT>
void loadSDPAFormatFile3(std::istream &is, LMII &lmi, VT &objectiveFunction) {
    std::string line;
    std::string::size_type sz;

    std::getline(is, line, '\n');

    //skip comments
    while (isCommentLine3(line)) {
        std::getline(is, line, '\n');
    }

    //read variables number
    int variablesNum = fetchNumber3(line);

    if (std::getline(is, line, '\n').eof())
        throw 1;

    //read number of blocks
    int blocksNum = fetchNumber3(line);

    if (std::getline(is, line, '\n').eof())
        throw 1;

    //read block structure vector
    listVector3 blockStructure = readVector3(line); //TODO different if we have one block

    if (blockStructure.size() != blocksNum)
        throw 1;

    if (std::getline(is, line, '\n').eof())
        throw 1;

    //read constant vector
    listVector3 constantVector = readVector3(line);

    while  (constantVector.size() < variablesNum) {
        if (std::getline(is, line, '\n').eof())
            throw 1;
        listVector3 t = readVector3(line);
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

                    listVector3 vec = readVector3(line);

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

                    listVector3 vec = readVector3(line);

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
Rcpp::NumericMatrix sample_points(Rcpp::Nullable<Rcpp::CharacterVector> file = R_NilValue,
                                  Rcpp::Nullable<std::string> distribution = R_NilValue,
                               Rcpp::Nullable<unsigned int> N = R_NilValue,
                               Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                                  Rcpp::Nullable<double> Temperature = R_NilValue,
                               Rcpp::Nullable<unsigned int> random_walk = R_NilValue){


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
    //std::string sdp = "sdp_prob"+bar+std::to_string(Rcpp::as<int>(nn))+bar+std::to_string(Rcpp::as<int>(mm))+txt;

    //std::cout<<"reading... "<<std::endl;

    inp.open(Rcpp::as<std::string>(file),std::ifstream::in);
    lmi Slmi;
    VT c;
    loadSDPAFormatFile3<MT>(inp, Slmi, c);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    

    spectaedro SP(Slmi);//, SP2;
    unsigned int n = SP.dimension();
    
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);
    
    Point p(n);
    std::list<Point> randPoints;
    NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0, diam, radius;
    bool round = false;
    NT T=1.0;
    int walkL = 1, NN=100;
    if(walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);
    if(N.isNotNull()) NN = Rcpp::as<unsigned int>(N);
    if(Temperature.isNotNull()) T = Rcpp::as<NT>(Temperature);

    SP.ComputeInnerBall(diam, radius);

    std::pair<Point,NT> InnerB;
    InnerB.first = p;

    vars<NT, RNGType> var(0,n, 1, 1,0.0,0.1,0,0.0,0, radius,diam,0,rng,urdist,urdist1,
                          -1.0,true,false,round,false,false,false,false,false, true);
    var.che_rad = radius;
    var.diameter = diam;

    if(!distribution.isNotNull() || Rcpp::as<std::string>(distribution).compare(std::string("uniform"))==0) {

        if (!random_walk.isNotNull() || Rcpp::as<std::string>(random_walk).compare(std::string("Billiard")) == 0) {


            spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
            settings.LMIatP = SP.getLMI().getA0();
            p = Point(n);

            rand_point_generator_spec(SP, p, NN, walkL, randPoints, var, settings);

        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("RDHR")) == 0) {
            p = Point(n);
            for (int j = 0; j < NN; ++j) {
                for (int k = 0; k < walkL; ++k) {
                    hit_and_run_spec(p, SP, var);
                }
                randPoints.push_back(p);
            }
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("CDHR")) == 0) {
            Point v(n);
            p = Point(n);
            int rand_coord;
            for (int j = 0; j < NN; ++j) {
                for (int k = 0; k < walkL; ++k) {
                    v = Point(n);
                    rand_coord = uidist(rng);
                    v.set_coord(rand_coord, 1.0);//(rand_coord) = 1.0;
                    std::pair <NT, NT> dbpair = SP.boundaryOracle(p.get_coefficients(), v.get_coefficients());
                    double min_plus = dbpair.first;
                    double max_minus = dbpair.second;
                    Point b1 = (min_plus * v) + p;
                    Point b2 = (max_minus * v) + p;
                    double lambda = urdist(rng);
                    p = (lambda * b1);
                    p = ((1 - lambda) * b2) + p;
                }
                randPoints.push_back(p);
            }
        }
    } else if (Rcpp::as<std::string>(distribution).compare(std::string("boltzmann"))==0){

        if (!random_walk.isNotNull() || Rcpp::as<std::string>(random_walk).compare(std::string("HMC")) == 0) {
            spectaedro::BoundaryOracleBoltzmannHMCSettings settings2;
            settings2.first = true;
            settings2.epsilon = 0.0001;
            Point cc(c);
            p = Point(n);
            for (int i = 0; i < NN; ++i) {
                for (int j = 0; j < walkL; ++j) {
                    HMC_boltzmann_reflections(SP, p, diam, var, cc, T, settings2);
                }
                randPoints.push_back(p);
            }
        } else if(Rcpp::as<std::string>(random_walk).compare(std::string("RDHR")) == 0){
            Point cc(c);
            p = Point(n);
            for (int i = 0; i < NN; ++i) {
                for (int j = 0; j < walkL; ++j) {
                    hit_and_run_Boltzmann_spec(p, SP, var, cc, T);
                }
                randPoints.push_back(p);
            }
        }else {
            throw Rcpp::exception("Wrong input!");
        }
    } else {
        throw Rcpp::exception("Wrong input!");
    }


    MT RetMat(n, NN);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = (*rpit).getCoefficients();

    return Rcpp::wrap(RetMat);


}