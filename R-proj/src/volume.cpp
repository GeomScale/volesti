// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iterator>
//#include <fstream>
#include <vector>
#include <list>
//#include <algorithm>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <random.hpp>
#include <random/uniform_int.hpp>
#include <random/normal_distribution.hpp>
#include <random/uniform_real_distribution.hpp>
#include "vars.h"
//#include "ellipsoids.h"
#include "ballintersectconvex.h"
#include "samplers.h"
#include "rounding.h"
#include "rotating.h"
#include "gaussian_annealing.h"
//#include "sample_only.h"
#include "spectrahedron.h"
//#include "volume.h"
#include "ball_vol_spec.h"
#include "sdp_generator.h"



typedef std::string::iterator string_it2;
typedef std::list<double> listVector2;


char consumeSymbol2(string_it2 &at, string_it2 &end) {
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


bool isCommentLine2(std::string &line) {
    string_it2 at = line.begin();
    string_it2 end = line.end();

    char c = consumeSymbol2(at, end);

    return c == '"' || c == '*';
}


int fetchNumber2(std::string &string) {
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
listVector2 readVector2(std::string &string) {
    std::stringstream stream(string);
    listVector2 vector;
    double value;

    while (stream >> value) {
        vector.push_back(value);
    }

    return vector;
}


template <typename MT, typename LMII, typename VT>
void loadSDPAFormatFile2(std::istream &is, LMII &lmi, VT &objectiveFunction) {
    std::string line;
    std::string::size_type sz;

    std::getline(is, line, '\n');

    //skip comments
    while (isCommentLine2(line)) {
        std::getline(is, line, '\n');
    }

    //read variables number
    int variablesNum = fetchNumber2(line);

    if (std::getline(is, line, '\n').eof())
        throw 1;

    //read number of blocks
    int blocksNum = fetchNumber2(line);

    if (std::getline(is, line, '\n').eof())
        throw 1;

    //read block structure vector
    listVector2 blockStructure = readVector2(line); //TODO different if we have one block

    if (blockStructure.size() != blocksNum)
        throw 1;

    if (std::getline(is, line, '\n').eof())
        throw 1;

    //read constant vector
    listVector2 constantVector = readVector2(line);

    while  (constantVector.size() < variablesNum) {
        if (std::getline(is, line, '\n').eof())
            throw 1;
        listVector2 t = readVector2(line);
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

                    listVector2 vec = readVector2(line);

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

                    listVector2 vec = readVector2(line);

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


//' The main function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope)
//'
//' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which correspond to the vertices of P. A zonotope is desrcibed by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param walk_step Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for SequenceOfBalls and \eqn{1} for CoolingGaussian.
//' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{1} for SequenceOfBalls and \eqn{0.1} for CoolingGaussian.
//' @param InnerBall Optional. A \eqn{d+1} vector that contains an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.
//' @param Algo Optional. A string that declares which algorithm to use: a) \code{'SoB'} for SequenceOfBalls or b) \code{'CG'} for CoolingGaussian.
//' @param WalkType Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run or c) \code{'BW'} for Ball Walk. The default walk is \code{'CDHR'}.
//' @param rounding Optional. A boolean parameter for rounding. The default value is \code{FALSE}.
//' @param Parameters Optional. A list for the parameters of the algorithms:
//' \itemize{
//' \item{\code{Window} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
//'  \item{\code{C} }{ A constant for the lower bound of \eqn{variance/mean^2} in schedule annealing of CG algorithm. The default value is \eqn{2}.}
//'  \item{\code{N} }{ The number of points we sample in each step of schedule annealing in CG algorithm. The default value is \eqn{500C + dimension^2 / 2}.}
//'  \item{\code{ratio} }{ Parameter of schedule annealing of CG algorithm, larger ratio means larger steps in schedule annealing. The default value is \eqn{1 - 1/dimension}.}
//'  \item{\code{frac} }{ The fraction of the total error to spend in the first gaussian in CG algorithm. The default value is \eqn{0.1}.}
//'  \item{\code{BW_rad} }{ The radius for the ball walk. The default value is \eqn{4r/dimension}, where \eqn{r} is the radius of the inscribed ball of the polytope.}
//' }
//'
//' @references \cite{I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2014.},
//' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
//'
//'
//' @return The approximation of the volume of a convex polytope.
//' @examples
//' # calling SOB algorithm for a H-polytope (2d unit simplex)
//' P = GenSimplex(2,'H')
//' vol = volume(P)
//'
//' # calling CG algorithm for a V-polytope (3d simplex)
//' P = GenSimplex(2,'V')
//' vol = volume(P, Algo = "CG")
//'
//' # calling CG algorithm for a 2-dimensional zonotope defined as the Minkowski sum of 4 segments
//' Z = GenZonotope(2, 4)
//' vol = volume(Z, WalkType = "RDHR", walk_step = 5)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector volume (Rcpp::Nullable<Rcpp::CharacterVector> file = R_NilValue,
               Rcpp::Nullable<unsigned int> walk_step = R_NilValue,
               Rcpp::Nullable<double> error = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> InnerBall = R_NilValue,
               Rcpp::Nullable<std::string> Algo = R_NilValue,
               Rcpp::Nullable<std::string> WalkType = R_NilValue,
               Rcpp::Nullable<bool> rounding = R_NilValue,
               Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue) {

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
    //std::string sdp = "sdp_prob"+bar+std::to_string(Rcpp::as<int>(d))+bar+std::to_string(Rcpp::as<int>(num))+txt;

    //std::cout<<"reading... "<<sdp<<std::endl;

    inp.open(Rcpp::as<std::string>(file),std::ifstream::in);
    lmi Slmi;
    VT c;
    loadSDPAFormatFile2<MT>(inp, Slmi, c);

    spectaedro SP(Slmi);//, SP2;
    unsigned int n = SP.dimension();

    bool CG, cdhr = true, rdhr = false, ball_walk = false, round=false, bref = false;
    unsigned int win_len = 4*n*n+500, N = 500 * 2 +  n * n / 2;

    double C = 2.0, ratio = 1.0-1.0/(NT(n)), frac = 0.1, e, delta = -1.0;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars_ban <NT> var_ban(0.1, 0.15, 0.75, 0.0, 0.2, 250, 125, 10, false);
    std::pair<Point,NT> InnerB;
    Point p(n);
    NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0;
    InnerB.first = p;// = SP.ComputeInnerBall(diam_spec);

    if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("bref")) {
        bref = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(Parameters)["bref"]);
    }

    round = (!rounding.isNotNull()) ? false : Rcpp::as<bool>(rounding);

    vars<NT, RNGType> var(0,n, 1, 1,0.0,0.1,0,0.0,0, InnerB.second,diam_spec,0,bref,rng,urdist,urdist1,
                          delta,true,false,round,false,false,false,false,false, true);
    spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
    settings.LMIatP = SP.getLMI().getA0();



    std::cout<<"\ninitializations ok.."<<std::endl;
    //vars<NT, RNGType> var2 = var;
    preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
    std::cout<<"preproccessing ok.."<<std::endl;
    InnerB.second = rad;

    var.steps=0;
    vol_spec = round_value * volesti_ball_ann(SP, var, var_ban, settings, InnerB, nballs2, false);

    Rcpp::NumericVector res(2);
    res[0] = vol_spec;
    res[1] = var.steps;

    return res;
}
