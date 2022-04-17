#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "chrono"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> spectrahedron;
typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;



std::vector<MT> Condition1(){
    std::vector<MT> list_Mat;
    //matrix
    /*
    -1,0,0
    0,-2,1
    0,1,-2
    */
    MT A;
    A = MT::Zero(3, 3);
    A(0,0)=-1;
    A(1,1)=-2;A(1,2)=1;
    A(2,1)=1;A(2,2)=-2;
    list_Mat.push_back(A);
    //matrix
    /*
    -1,0,0
    0,0,1
    0,1,0
    */
    A = MT::Zero(3, 3);
    A(0,0)=-1;
    A(1,2)=1;
    A(2,1)=1;
    list_Mat.push_back(A);
    //matrix
    /*
    0,0,-1
    0,0,0
    -1,0,0
    */
    A = MT::Zero(3, 3);
    A(0,2)=-1;
    A(2,0)=-1;
    list_Mat.push_back(A);


    return list_Mat;
}
std::vector<MT> Condition2(){
    std::vector<MT> list_Mat;
    //matrix
    /*
    -5,0
    0,-5
    */
    MT A;
    A = MT::Zero(2, 2);
    A(0,0)=-5;
    A(1,1)=-5;
    list_Mat.push_back(A);
    //matrix
    /*
    1,0
    0,-1
    */
    A = MT::Zero(2, 2);
    A(0,0)=1;
    A(1,1)=-1;
    list_Mat.push_back(A);
    //matrix
    /*
    0,0,-1
    0,0,0
    -1,0,0
    */
    A = MT::Zero(2, 2);
    A(0,1)=4;
    A(1,0)=4;
    list_Mat.push_back(A);


    return list_Mat;
}


Point getDirection(unsigned int const& dim, RNGType &rng, bool normalize=true){
    double normal = 0.;
    Point p(dim);
    double* data = p.pointerToData();

    for (unsigned int i=0; i<dim; ++i){
        *data = rng.sample_ndist();
        normal += *data * *data;
        data++;
    }

    normal = 1./std::sqrt(normal);
    if (normalize) p *= normal;
    return p;
}

/*Pick a uniform direction on the halfplane of the normal vector */
void unifDirection(spectrahedron P, Point &p, Point &v, RNGType &rng){
    unsigned int d = P.dimension();

    VT nt(d);
    (P.lmi).normalizedDeterminantGradient(p.getCoefficients(), (P.precomputedValues).eigenvector, nt);
    v=getDirection(d,rng);

    if(v.dot(nt)>0){v=-1*v;}
}
/*Stochastic Billiards algorithm */
std::vector<Point> StochasticBilliard(spectrahedron &P, Point& q, unsigned int num_points, RNGType &rng){
    unsigned int d = P.dimension();
    double  tau;
    Point p = q, v;
    std::vector<Point> randPoints;
        v = getDirection(d, rng);
        int it = 0;
        //For each sample
        while (it < num_points)
        {
          //Find the intersection
            tau=P.positiveLinearIntersection(p.getCoefficients(), v.getCoefficients());
            randPoints.push_back(p+tau*v);
            tau = 0.995 * tau;
            //Move to it
            p += tau * v;
            //pick a uniform direction
            unifDirection(P, p, v,rng);
            it++;
        }
    return randPoints;
}
//Entry function
std::vector<Point> boundarySampler(unsigned int d, unsigned int num_points){
    RNGType rng(d);
    Point p(d);
    // Create a spectrahedron
    std::vector<MT> lmi_mat = Condition1();
    LMI<double, MT, VT> lmi(lmi_mat);
    spectrahedron spectra(lmi);
    spectra.set_interior_point(p);
    return StochasticBilliard(spectra, p, num_points, rng);

}

int main(int argc, char** argv){
    if(argc!=2){std::cout<<"usage: ./sampler numberOfSamples";}
    srand((unsigned) time(NULL));
    int d = 2;
    int numberOfSamples=atoi(argv[1]);
    std::vector<Point> points = boundarySampler(d, numberOfSamples);
    for(int i=0;i<numberOfSamples;i++){
      for(int j=0;j<d;j++){
        std::cout<<(points[i])[j]<<" ";
      }
      std::cout<<"\n";
    }
    return 0;
}
