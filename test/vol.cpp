// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#include "Eigen/Eigen"
#define VOLESTI_DEBUG
#include <fstream>
#include <iterator>
//#include <fstream>
#include <vector>
#include <list>
//#include <algorithm>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
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
//#include "sdp_generator.h"
//#include "sample_only.h"
//#include "exact_vols.h"

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

template <typename FT>
FT factorial(FT n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


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

// Approximating the volume of a convex polytope or body 
// can also be used for integration of concave functions.
// The user should provide the appropriate membership 
// oracles.

int main(const int argc, const char** argv)
{
	//Deafault values
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT, NT, VT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <MT, VT> lmi;
    typedef Spectrahedron <lmi, Point> spectaedro;

    int n, nexp=1, n_threads=1, W;
    int walk_len=1,N=100,max_iter=20, nsam = 100;
    NT e=1;
    NT exactvol(-1.0), a=0.5;
    VT c;
    bool verbose=false, 
	 rand_only=false, 
	 round_only=false,
	 file=false, 
	 round=false, 
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         ball_walk=false,
         ball_rad=false,
         experiments=true,
         annealing = false,
         Vpoly=false,
         Zono=false,
         cdhr=false,
         rdhr=false,
         exact_zono = false,
         billiard=true,
         hmc = true,
            sample_only = false,
            sdp=false,
            boltz = false,
         gaussian_sam = false;
    spectaedro SP;

    //this is our polytope
    //Hpolytope HP;
    //Vpolytope VP; // RNGType only needed for the construction of the inner ball which needs randomization
    //Zonotope  ZP;

    // parameters of CV algorithm
    bool user_W=false, user_N=false, user_ratio=false;
    NT ball_radius=0.0, T = 1.0;
    NT C=2.0,ratio,frac=0.1,delta=-1.0,error=0.2;

  if(argc<2){
    std::cout<<"Use -h for help"<<std::endl;
    exit(-2);
  }
  
  //parse command line input vars
  for(int i=1;i<argc;++i){
      bool correct = false;





      if(!strcmp(argv[i],"-sample")){
          sample_only = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-hmc")){
          hmc = true;
          billiard = false;
          correct = true;
      }
      if(!strcmp(argv[i],"-rdhr")){
          rdhr = true;
          hmc=false;
          billiard = false;
          correct = true;
      }
      if(!strcmp(argv[i],"-cdhr")){
          cdhr = true;
          billiard = false;
          correct = true;
      }
      if(!strcmp(argv[i],"-billiard")){
          billiard = true;
          correct = true;
      }


      //if(!strcmp(argv[i],"-max_iter")){
          //max_iter = atof(argv[++i]);
          //correct = true;
      //}
      if(!strcmp(argv[i],"-N")){
          N = atof(argv[++i]);
          max_iter = N;
          user_N = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-temperature")){
          T = atof(argv[++i]);
          correct = true;
      }

      if(!strcmp(argv[i],"-sdp")){
          sdp = true;
          correct = true;
      }
      if(!strcmp(argv[i],"-boltz")){
          boltz = true;
          correct = true;
      }
      //reading from file
      if(!strcmp(argv[i],"-file")){
          std::ifstream inp;


          //std::cout<<"reading spactrahedra... "<<std::endl;

          inp.open(argv[++i],std::ifstream::in);
          lmi Slmi;

          loadSDPAFormatFile3<MT>(inp, Slmi, c);
          spectaedro SP2(Slmi);//, SP2;
          n = SP2.dimension();
          SP = SP2;
          correct = true;
      }

      //reading linear extensions and order polytopes

      if(!strcmp(argv[i],"-w")||!strcmp(argv[i],"-walk_length")){
          walk_len = atof(argv[++i]);
          user_walk_len = true;
          correct = true;
      }

      if(correct==false){
          std::cerr<<"unknown parameters \'"<<argv[i]<<
                     "\', try "<<argv[0]<<" --help"<<std::endl;
          exit(-2);
      }
      
  }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

  if(!sample_only && !sdp) {
      //cdhr = true; rdhr = false; ball_walk = false; round=false;
      //int win_len = 4*n*n+500;
      //N = 500 * 2 +  n * n / 2;

      //C = 2.0, ratio = 1.0-1.0/(NT(n)), frac = 0.1, e, delta = -1.0;

      vars_ban <NT> var_ban(0.05, 0.1, 0.75, 0.0, 0.2, 150, 125, 10, false);
      std::pair<Point,NT> InnerB;
      Point p(n);
      NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0;
      InnerB.first = p;// = SP.ComputeInnerBall(diam_spec);

      vars<NT, RNGType> var(0,n, 1, 1,0.0,0.1,0,0.0,0, InnerB.second,diam_spec,rng,urdist,urdist1,
                            delta,true,false,round,false,false,false,false,false, true);
      spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
      settings.LMIatP = SP.getLMI().getA0();

      //std::cout<<"\ninitializations ok.."<<std::endl;
      preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
      //std::cout<<"preproccessing ok.."<<std::endl;
      InnerB.second = rad;

      vol_spec = volesti_ball_ann(SP, var, var_ban, settings, InnerB, nballs2, false);

      std::cout<<"volume = "<<vol_spec<<std::endl;
      return -1.0;
  } else if(sample_only) {
      Point p(n);
      std::list<Point> randPoints;
      NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0, diam, radius;
      bool round = false;

      SP.ComputeInnerBall(diam, radius);

      std::pair<Point,NT> InnerB;
      InnerB.first = p;

      vars<NT, RNGType> var(0,n, 1, 1,0.0,0.1,0,0.0,0, radius,diam,rng,urdist,urdist1,
                            -1.0,true,false,round,false,false,false,false,false, true);
      var.che_rad = radius;
      var.diameter = diam;

      if(!boltz) {
          if (billiard) {


              spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
              settings.LMIatP = SP.getLMI().getA0();
              p = Point(n);

              rand_point_generator_spec(SP, p, N, walk_len, randPoints, var, settings);

          } else if (hmc) {
              spectaedro::BoundaryOracleBoltzmannHMCSettings settings2;
              settings2.first = true;
              settings2.epsilon = 0.0001;
              Point cc(c);
              p = Point(n);
              for (int i = 0; i < N; ++i) {
                  for (int j = 0; j < walk_len; ++j) {
                      HMC_boltzmann_reflections(SP, p, diam, var, cc, T, settings2);
                  }
                  randPoints.push_back(p);
              }
          } else if (rdhr) {
              p = Point(n);
              for (int j = 0; j < N; ++j) {
                  for (int k = 0; k < walk_len; ++k) {
                      hit_and_run_spec(p, SP, var);
                  }
                  randPoints.push_back(p);
              }
          } else if (cdhr) {
              Point v(n);
              p = Point(n);
              int rand_coord;
              for (int j = 0; j < NN; ++j) {
                  for (int k = 0; k < walk_len; ++k) {
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
      } else {
          if (hmc) {
              spectaedro::BoundaryOracleBoltzmannHMCSettings settings2;
              settings2.first = true;
              settings2.epsilon = 0.0001;
              Point cc(c);
              p = Point(n);
              for (int i = 0; i < N; ++i) {
                  for (int j = 0; j < walk_len; ++j) {
                      HMC_boltzmann_reflections(SP, p, diam, var, cc, T, settings2);
                  }
                  randPoints.push_back(p);
              }
          } else {
              Point cc(c);
              p = Point(n);
              for (int i = 0; i < N; ++i) {
                  for (int j = 0; j < walk_len; ++j) {
                      hit_and_run_Boltzmann_spec(p, SP, var, cc, T);
                  }
                  randPoints.push_back(p);
              }
          }
      }

      for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++) {
          (*rpit).print();
      }

  } else {
      bool round = false;
      //std::pair <Point, NT> InnerB;
      Point p(n);
      NT nballs2, diam, vol_spec, radius, round_value = 1.0;
      SP.ComputeInnerBall(diam, radius);

      std::pair<Point,NT> InnerB;
      InnerB.first = p;

      vars<NT, RNGType> var(0,n, 1, 1,0.0,0.1,0,0.0,0, radius,diam,rng,urdist,urdist1,
                            -1.0,true,false,round,false,false,false,false,false, true);
      var.che_rad = radius;
      var.diameter = diam;

      //std::list <Point> randPoints, randPoints2;
      //spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
      //settings.LMIatP = SP.getLMI().getA0();
      //preproccess_spectrahedron(SP, p, var, settings, round_value, diam, rad, round);
      //settings.LMIatP = SP.getLMI().getA0();
      //p = Point(n);

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
      std::vector<NT>retvec;


      //if(max_iterations.isNotNull()) max_iter = Rcpp::as<unsigned int>(max_iterations);
      //Rcpp::NumericVector retvec(NN);

      if(hmc){


          //std::cout << "HMC"<< std::endl;
          //std::cout << optimal_solution<< std::endl;
          //s/td::cout << err<< std::endl;
          //std::cout << (std::abs(min_val - optimal_solution) / std::abs(optimal_solution)) << std::endl;

          for (int i = 0; i < max_iter; ++i) {
              for (int j = 0; j < walk_len; ++j) {
                  HMC_boltzmann_reflections(SP, p, diam, var, cc, T, settings2);
              }
              retvec.push_back(c.dot(p.getCoefficients()));
              T = T * tred;

          }
      } else if(rdhr) {


          //std::cout << "Hit and Run" << std::endl;

          for (int i = 0; i < max_iter; ++i) {
              for (int j = 0; j < walk_len; ++j) {
                  hit_and_run_Boltzmann_spec(p, SP, var, cc, T);
              }
              retvec.push_back(c.dot(p.getCoefficients()));
              T = T * tred;

          }
      }

      std::cout<<"The values of the objective function in each iteration:\n"<<std::endl;
      //std::cout<<", ";
      for (int k = 0; k < max_iter; ++k) {
          if (k>0) std::cout<<", ";
          std::cout<<float(retvec[k]);
      }
      std::cout<<"\n";

  }


  return 0;
}
