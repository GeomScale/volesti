// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "sampling/sampling.hpp"

template <typename Polytope, typename RNGType, typename PointList, typename NT, typename Point>
void sample_from_polytope(Polytope &P, RNGType &rng, PointList &randPoints, unsigned int const& walkL, unsigned int const& numpoints,
        bool const& gaussian, NT const& a, NT const& L, bool const& boundary, Point const& StartingPoint, unsigned int const& nburns,
        bool const& set_L, bool const& cdhr, bool const& rdhr, bool const& billiard, bool const& ball_walk)
{
    if (boundary) {
        if (cdhr) {
            uniform_sampling_boundary <BCDHRWalk>(randPoints, P, rng, walkL, numpoints,
                     StartingPoint, nburns);
        } else {
            uniform_sampling_boundary <BRDHRWalk>(randPoints, P, rng, walkL, numpoints,
                     StartingPoint, nburns);
        }
    } else if (cdhr) {
        if (gaussian) {
            gaussian_sampling<GaussianCDHRWalk>(randPoints, P, rng, walkL, numpoints,
                                             a, StartingPoint, nburns);
        } else {
            uniform_sampling<CDHRWalk>(randPoints, P, rng, walkL, numpoints,
                                             StartingPoint, nburns);
        }
    } else if (rdhr){
        if (gaussian) {
            gaussian_sampling<GaussianRDHRWalk>(randPoints, P, rng, walkL, numpoints,
                                             a, StartingPoint, nburns);
        } else {
            uniform_sampling<RDHRWalk>(randPoints, P, rng, walkL, numpoints,
                                             StartingPoint, nburns);
        }
    } else if (billiard) {
        if (set_L) {
            BilliardWalk WalkType(L);
            uniform_sampling(randPoints, P, rng, WalkType, walkL, numpoints, StartingPoint, nburns);
        } else {
            uniform_sampling<BilliardWalk>(randPoints, P, rng, walkL, numpoints,
                     StartingPoint, nburns);
        }
    } else {
        if (set_L) {
            if (gaussian) {
                GaussianBallWalk WalkType(L);
                gaussian_sampling(randPoints, P, rng, WalkType, walkL, numpoints, a,
                                       StartingPoint, nburns);
            } else {
                BallWalk WalkType(L);
                uniform_sampling(randPoints, P, rng, WalkType, walkL, numpoints,
                                       StartingPoint, nburns);
            }
        } else {
            if (gaussian) {
                gaussian_sampling<GaussianBallWalk>(randPoints, P, rng, walkL, numpoints,
                                                 a, StartingPoint, nburns);
            } else {
                uniform_sampling<BallWalk>(randPoints, P, rng, walkL, numpoints,
                                                 StartingPoint, nburns);
            }
        }
    }
}


//' Sample uniformly or normally distributed points from a convex Polytope (H-polytope, V-polytope, zonotope or intersection of two V-polytopes).
//'
//' Sample n points with uniform or multidimensional spherical gaussian -with a mode at any point- as the target distribution.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope or (d) VpolytopeIntersection.
//' @param n The number of points that the function is going to sample from the convex polytope.
//' @param random_walk Optional. A list that declares the random walk and some related parameters as follows:
//' \itemize{
//' \item{\code{walk} }{ A string to declare the random walk: i) \code{'CDHR'} for Coordinate Directions Hit-and-Run, ii) \code{'RDHR'} for Random Directions Hit-and-Run, iii) \code{'BaW'} for Ball Walk, iv) \code{'BiW'} for Billiard walk, v) \code{'BCDHR'} boundary sampling by keeping the extreme points of CDHR or vi) \code{'BRDHR'} boundary sampling by keeping the extreme points of RDHR. The default walk is \code{'BiW'} for the uniform distribution or \code{'CDHR'} for the Gaussian distribution.}
//' \item{\code{walk_length} }{ The number of the steps per generated point for the random walk. The default value is 1.}
//' \item{\code{nburns} }{ The number of points to burn before start sampling.}
//' \item{\code{starting_point} }{ A \eqn{d}-dimensional numerical vector that declares a starting point in the interior of the polytope for the random walk. The default choice is the center of the ball as that one computed by the function \code{inner_ball()}.}
//' \item{\code{BaW_rad} }{ The radius for the ball walk.}
//' \item{\code{L} }{ The maximum length of the billiard trajectory.}
//' \item{\code{seed} }{ A fixed seed for the number generator.}
//' }
//' @param distribution Optional. A list that declares the target density and some related parameters as follows:
//' \itemize{
//' \item{\code{density} }{ A string: (a) \code{'uniform'} for the uniform distribution or b) \code{'gaussian'} for the multidimensional spherical distribution. The default target distribution is uniform.}
//' \item{\code{variance} }{ The variance of the multidimensional spherical gaussian. The default value is 1.}
//'  \item{\code{mode} }{ A \eqn{d}-dimensional numerical vector that declares the mode of the Gaussian distribution. The default choice is the center of the as that one computed by the function \code{inner_ball()}.}
//' }
//'
//' @return A \eqn{d\times n} matrix that contains, column-wise, the sampled points from the convex polytope P.
//' @examples
//' # uniform distribution from the 3d unit cube in H-representation using ball walk
//' P = gen_cube(3, 'H')
//' points = sample_points(P, n = 100, random_walk = list("walk" = "BaW", "walk_length" = 5))
//'
//' # gaussian distribution from the 2d unit simplex in H-representation with variance = 2
//' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
//' b = c(0,0,1)
//' P = Hpolytope(A = A, b = b)
//' points = sample_points(P, n = 100, distribution = list("density" = "gaussian", "variance" = 2))
//'
//' # uniform points from the boundary of a 2-dimensional random H-polytope
//' P = gen_rand_hpoly(2,20)
//' points = sample_points(P, n = 100, random_walk = list("walk" = "BRDHR"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_points(Rcpp::Reference P,
                                  int n,
                                  Rcpp::Nullable<Rcpp::List> random_walk = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> distribution = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point    Point;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope, RNGType> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    unsigned int dim, walkL = 1, type_num;

    std::string type = Rcpp::as<std::string>(P.slot("type"));

    if (type.compare(std::string("Hpolytope")) == 0) {
        dim = Rcpp::as<MT>(P.slot("A")).cols();
        type_num = 1;
    } else if (type.compare(std::string("Vpolytope")) == 0) {
        dim = Rcpp::as<MT>(P.slot("V")).cols();
        type_num = 2;
    } else if (type.compare(std::string("Zonotope")) == 0) {
        dim = Rcpp::as<MT>(P.slot("G")).cols();
        type_num = 3;
    } else if (type.compare(std::string("VpolytopeIntersection")) == 0) {
        dim = Rcpp::as<MT>(P.slot("V1")).cols();
        type_num = 4;
    } else {
        throw Rcpp::exception("Unknown polytope representation!");
    }

    RNGType rng(dim);
    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("seed")) {
        unsigned seed2 = Rcpp::as<double>(Rcpp::as<Rcpp::List>(random_walk)["seed"]);
        rng.set_seed(seed2);
    }

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;
    InterVP VPcVP;

    unsigned int numpoints, nburns = 0;
    NT radius = 1.0, L;
    bool set_mode = false, cdhr = false, rdhr = false, ball_walk = false, gaussian = false,
          billiard = false, boundary = false, set_starting_point = false, set_L = false;
    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall;
    Point mode(dim);

    numpoints = n;
    if (numpoints <= 0) throw Rcpp::exception("The number of samples has to be a positice integer!");

    if (!distribution.isNotNull() || !Rcpp::as<Rcpp::List>(distribution).containsElementNamed("density")) {
        billiard = true;
    } else if (
            Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(distribution)["density"]).compare(std::string("uniform")) == 0) {
        billiard = true;
    } else if (
            Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(distribution)["density"]).compare(std::string("gaussian")) == 0) {
        gaussian = true;
    } else {
        throw Rcpp::exception("Wrong distribution!");
    }

    if (Rcpp::as<Rcpp::List>(distribution).containsElementNamed("mode")) {
        if (!gaussian) throw Rcpp::exception("Mode is given only for Gaussian sampling!");
        if (Rcpp::as<VT>(Rcpp::as<Rcpp::List>(distribution)["mode"]).size() != dim) {
            throw Rcpp::exception("Mode has to be a point in the same dimension as the polytope P");
        } else {
            set_mode = true;
            VT temp = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(distribution)["mode"]);
            mode = Point(dim, std::vector<NT>(&temp[0], temp.data() + temp.cols() * temp.rows()));
        }
    }

    NT a = 0.5;
    if (Rcpp::as<Rcpp::List>(distribution).containsElementNamed("variance")) {
        a = 1.0 / (2.0 * Rcpp::as<NT>(Rcpp::as<Rcpp::List>(distribution)["variance"]));
        if (!gaussian) {
            Rcpp::warning("The variance can be set only for Gaussian sampling!");
        } else if (a <= 0.0) {
            throw Rcpp::exception("The variance has to be positive!");
        }
    }

    if (!random_walk.isNotNull() || !Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("walk")) {
        if (gaussian) {
            if (type_num == 1) {
                cdhr = true;
            } else {
                rdhr = true;
            }
        } else {
            billiard = true;
        }
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("CDHR")) == 0) {
        cdhr = true;
        billiard = false;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("RDHR")) == 0) {
        rdhr = true;
        billiard = false;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BaW")) == 0) {
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("BaW_rad")) {
            L = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(random_walk)["BaW_rad"]);
            set_L = true;
            if (L<=0.0) throw Rcpp::exception("BaW diameter must be a postitive number!");
        }
        ball_walk = true;
        billiard = false;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BiW")) == 0) {
        if (gaussian) throw Rcpp::exception("Billiard walk can be used only for uniform sampling!");
        billiard = true;
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("L")) {
            L = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(random_walk)["L"]);
            set_L = true;
            if (L<=0.0) throw Rcpp::exception("L must be a postitive number!");
        }
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BRDHR")) == 0) {
        if (gaussian) throw Rcpp::exception("Gaussian sampling from the boundary is not supported!");
        rdhr = true;
        boundary = true;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BCDHR")) == 0) {
        if (gaussian) throw Rcpp::exception("Gaussian sampling from the boundary is not supported!");
        cdhr = true;
        boundary = true;
    } else {
        throw Rcpp::exception("Unknown walk type!");
    }

    Point StartingPoint;
    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("starting_point")) {
        if (Rcpp::as<VT>(Rcpp::as<Rcpp::List>(random_walk)["starting_point"]).size() != dim) {
            throw Rcpp::exception("Starting Point has to lie in the same dimension as the polytope P");
        } else {
            set_starting_point = true;
            VT temp = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(random_walk)["starting_point"]);
            StartingPoint = Point(dim, std::vector<NT>(&temp[0], temp.data() + temp.cols() * temp.rows()));
        }
    }

    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("walk_length")) {
        walkL = Rcpp::as<int>(Rcpp::as<Rcpp::List>(random_walk)["walk_length"]);
        if (walkL <= 0) {
            throw Rcpp::exception("The walk length has to be a positive integer!");
        }
    }


    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("nburns")) {
        nburns = Rcpp::as<int>(Rcpp::as<Rcpp::List>(random_walk)["nburns"]);
        if (nburns < 0) {
            throw Rcpp::exception("The number of points to burn before sampling has to be a positive integer!");
        }
    }

    switch(type_num) {
        case 1: {
            // Hpolytope
            HP.init(dim, Rcpp::as<MT>(P.slot("A")),
                    Rcpp::as<VT>(P.slot("b")));

            if (!set_starting_point || (!set_mode && gaussian)) {
                InnerBall = HP.ComputeInnerBall();
                if (!set_starting_point) StartingPoint = InnerBall.first;
                if (!set_mode && gaussian) mode = InnerBall.first;
            }
            if (HP.is_in(StartingPoint) == 0) {
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            }
            HP.normalize();
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                HP.shift(mode.getCoefficients());
            }
            break;
        }
        case 2: {
            // Vpolytope
            VP.init(dim, Rcpp::as<MT>(P.slot("V")),
                    VT::Ones(Rcpp::as<MT>(P.slot("V")).rows()));

            if (!set_starting_point || (!set_mode && gaussian)) {
                InnerBall = VP.ComputeInnerBall();
                if (!set_starting_point) StartingPoint = InnerBall.first;
                if (!set_mode && gaussian) mode = InnerBall.first;
            }
            if (VP.is_in(StartingPoint) == 0)
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                VP.shift(mode.getCoefficients());
            }
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(dim, Rcpp::as<MT>(P.slot("G")),
                    VT::Ones(Rcpp::as<MT>(P.slot("G")).rows()));

            if (!set_starting_point || (!set_mode && gaussian)) {
                InnerBall = ZP.ComputeInnerBall();
                if (!set_starting_point) StartingPoint = InnerBall.first;
                if (!set_mode && gaussian) mode = InnerBall.first;
            }
            if (ZP.is_in(StartingPoint) == 0)
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                ZP.shift(mode.getCoefficients());
            }
            break;
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1;
            Vpolytope VP2;
            VP1.init(dim, Rcpp::as<MT>(P.slot("V1")),
                     VT::Ones(Rcpp::as<MT>(P.slot("V1")).rows()));
            VP2.init(dim, Rcpp::as<MT>(P.slot("V2")),
                     VT::Ones(Rcpp::as<MT>(P.slot("V2")).rows()));
            VPcVP.init(VP1, VP2);

            if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
            InnerBall = VPcVP.ComputeInnerBall();
            if (!set_starting_point) StartingPoint = InnerBall.first;
            if (!set_mode && gaussian) mode = InnerBall.first;
            if (VPcVP.is_in(StartingPoint) == 0)
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                VPcVP.shift(mode.getCoefficients());
            }
            break;
        }
    }

    switch (type_num) {
        case 1: {
            sample_from_polytope(HP, rng, randPoints, walkL, numpoints, gaussian, a, L, boundary, StartingPoint, nburns,
                                 set_L, cdhr, rdhr, billiard, ball_walk);
            break;
        }
        case 2: {
            sample_from_polytope(VP, rng, randPoints, walkL, numpoints, gaussian, a, L, boundary, StartingPoint, nburns,
                                 set_L, cdhr, rdhr, billiard, ball_walk);
            break;
        }
        case 3: {
            sample_from_polytope(ZP, rng, randPoints, walkL, numpoints, gaussian, a, L, boundary, StartingPoint, nburns,
                                 set_L, cdhr, rdhr, billiard, ball_walk);
            break;
        }
        case 4: {
            sample_from_polytope(VPcVP, rng, randPoints, walkL, numpoints, gaussian, a, L, boundary, StartingPoint, nburns,
                                 set_L, cdhr, rdhr, billiard, ball_walk);
            break;
        }
    }

    if (numpoints % 2 == 1 && boundary) numpoints--;
    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
        if (gaussian) {

            RetMat.col(jj) = rpit->getCoefficients() + mode.getCoefficients();

        } else {
            RetMat.col(jj) = (*rpit).getCoefficients();
        }
    }

    return Rcpp::wrap(RetMat);

}
