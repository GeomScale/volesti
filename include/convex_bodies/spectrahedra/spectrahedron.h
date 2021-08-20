// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "LMI.h"
#include "chrono"


template <typename MT>
struct PrecomputedValuesLineIntersection {
    
};


template <typename NT>
struct PrecomputedValuesLineIntersection<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>> {

    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// These flags indicate whether the corresponding matrices are computed
    /// if yes, we can use them and not compute them fro scratch
    bool computed_B = false;
    bool computed_C = false;

    MT B, C;
    VT eigenvector;

    /// Sets all flags to false
    void resetFlags() {
        computed_C = computed_B = false;
    }

    void set_mat_size(int const& m) 
    {
        B.setZero(m, m);
        C.setZero(m, m);
        eigenvector.setZero(m);
    }
};


/// This class manipulates a spectrahedron, described by a LMI
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class Spectrahedron {
public:

    /// The numeric/matrix/vector types we use
    typedef NT NUMERIC_TYPE;
    typedef MT MATRIX_TYPE;
    typedef VT VECTOR_TYPE;

    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point PointType;

    /// The type of a pair of NT
    typedef std::pair<NT, NT> pairNT;

    /// Among successive calls of this class methods, we may need to pass data
    /// from one call to the next, to avoid repeating computations, or to efficiently update values
    /// Warning: this struct assists in many methods; perhaps for different methods use different instances
    struct PrecomputedValues {

        /// These flags indicate whether the corresponding matrices are computed
        /// if yes, we can use them and not compute them fro scratch
        bool computed_A = false;
        bool computed_C = false;
        bool computed_XY = false;

        /// The matrices the method positiveIntersection receives from its previous call
        /// if the flag first_positive_intersection is true.
        /// Matrix A is also used in coordinateIntersection
        MT A, B, C, X, Y;

        /// In method positive_intersect, the distance we are computing corresponds
        /// to the minimum positive eigenvalue of a quadratic eigenvalue problem.
        /// This will hold the eigenvector for that eigenvalue
        VT eigenvector;

        /// Sets all flags to false
        void resetFlags() {
            computed_XY = computed_C = computed_A = false;
        }
    };

    typedef PrecomputedValuesLineIntersection<MT> PrecomputedValuesLine;
    PrecomputedValuesLine precomputedVals;


    /// The dimension of the spectrahedron
    unsigned int d;
    VT grad;

    /// The linear matrix inequality that describes the spectrahedron
    LMI<NT, MT, VT> lmi;

    Spectrahedron() {}

    /// Creates a spectrahedron
    /// \param[in] lmi The linear matrix inequality that describes the spectrahedron
    Spectrahedron(const LMI<NT, MT, VT>& lmi) : lmi(lmi) {
        d = lmi.dimension();
        grad.setZero(d);
        precomputedVals.resetFlags();
        precomputedVals.set_mat_size(lmi.sizeOfMatrices());
    }


    /// Construct the quadratic eigenvalue problem \[At^2 + Bt + C \] for positive_intersect.
    /// A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// \param[in] a Input vector
    /// \param[in] b Input vector
    /// \param[in] c Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPositiveIntersection(const VT& a, const VT& b, const VT& c,
                                               PrecomputedValues& precomputedValues) {
        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_A) {
            lmi.evaluateWithoutA0(a, precomputedValues.A);
        }

        if (!precomputedValues.computed_C) {
            lmi.evaluate(c, precomputedValues.C);
        }

        // compute Matrix B
        lmi.evaluateWithoutA0(b, precomputedValues.B);
    }


     void createMatricesForPositiveIntersection(const VT& p, const VT& v) {
        
        // check if matrices B, C are ready if not compute them
        if (!precomputedValues.computed_C) 
        {
            lmi.evaluate(p, precomputedValues.C);
        }

        lmi.evaluateWithoutA0(v, precomputedValues.B);
    }

    /// Computes the distance d we must travel on the parametrized polynomial curve \[at^2 + bt + c \],
    /// assuming we start at t=0, and we start increasing t.
    /// We construct the quadratic eigenvalue problem \[At^2 + Bt + C \],
    /// where A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// Then we do a linearization and transform it to the generalized problem X+lY,
    /// which we pass to an external library.
    /// \param[in] a Input vector, the coefficient of t \[t^2\]
    /// \param[in] b Input vector, the coefficient of t
    /// \param[in] c Input Vector, the constant term
    /// \param[in, out] precomputedValues Data we move between successive method calls
    /// \returns The distance d
    NT positiveIntersection(VT const & a, VT const & b, VT const & c, PrecomputedValues& precomputedValues) 
    {
        unsigned int matrixDim = lmi.sizeOfMatrices();

        // create matrices A, B, C
        createMatricesForPositiveIntersection(a, b, c, precomputedValues);

        // get the minimum positive eigenvalue of At^2 + Bt + C
        EigenvaluesProblems<NT, MT, VT> quadraticEigenvaluesProblem;
        NT distance = quadraticEigenvaluesProblem.minPosQuadraticEigenvalue(precomputedValues.A, precomputedValues.B,
                                                                            precomputedValues.C, precomputedValues.X,
                                                                            precomputedValues.Y,
                                                                            precomputedValues.eigenvector,
                                                                            precomputedValues.computed_XY);
        return distance;
    }


    NT positiveIntersection(VT const & p, VT const & v) 
    {
        createMatricesForPositiveIntersection(p, v);
        EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;
        NT distance = EigenvaluesProblem.minPosLinearEigenvalue(precomputedVals.C, precomputedVals.B,
                                                                precomputedVals.eigenvector);
        return distance;
    }


    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType& r,
                                               PointType& v,
                                               VT&,
                                               VT& ,
                                               NT&,
                                               update_parameters&) 
    {
        NT pos_inter = positiveIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType& r,
                                               PointType& v,
                                               VT& ,
                                               VT& ,
                                               NT &,
                                               MT& ,
                                               update_parameters& ) 
    {
        NT pos_inter = positiveIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType& r,
                                                     PointType& v,
                                                     VT& ,
                                                     VT& ,
                                                     update_parameters&)
    {
        NT pos_inter = positiveIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    /// Computes the distance d one must travel on the line a + tb,
    /// assuming we start at t=0 and that b has zero everywhere and 1 in its i-th coordinate.
    /// We must solve the generalized eigenvalue problem A+tB, where A = lmi(a) and B=(lmi) - A0 = A_i
    /// If the flag precomputedValues,computed_A is true, the matrix A is not computed.
    /// \param[in] a Input vector
    /// \param[in] coordinate Indicator of the i-th coordinate, 1 <= coordinate <= dimension
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT coordinateIntersection(VT const & a, int const coordinate, PrecomputedValues& precomputedValues) {

        // prepare the generalized eigenvalue problem A+lB
        // we may not have to compute A!
        if (!precomputedValues.computed_A)
            lmi.evaluate(a, precomputedValues.A);

        EigenvaluesProblems<NT, MT, VT> eigenvaluesProblems;
        return eigenvaluesProblems.symGeneralizedProblem(precomputedValues.A, *(lmi.getMatrix(coordinate)));
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] point A point on the boundary of the spectrahedron
    /// \param[in] incomingDirection The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    /// \param[in] precomputedValues Must contain an eigenvalue needed to compute the reflection
    void computeReflection(VT const & point, VT const & incomingDirection, VT& reflectedDirection,
            PrecomputedValues& precomputedValues) {

        // get the gradient of the determinant of the lmi at point
        VT grad;
        lmi.normalizedDeterminantGradient(point, precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * incomingDirection.dot(grad);
        reflectedDirection = incomingDirection - dot * grad;
    }

    template <typename update_parameters>
    void compute_reflection(PointType &v, const PointType &r, update_parameters& ) {

        lmi.normalizedDeterminantGradient(r.getCoefficients(), precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * v.dot(grad);
        v += -dot * PointType(grad);
    }

    void compute_reflection(PointType &v, const PointType &r, int&) {

        lmi.normalizedDeterminantGradient(r.getCoefficients(), precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * v.dot(grad);
        v += -dot * PointType(grad);
    }


    /// \return The dimension of the spectrahedron
    unsigned int dimension() const {
        return d;
    }

    /// \return The LMI describing this spectrahedron
    LMI<NT, MT, VT> getLMI() const {
        return lmi;
    }

    /// Estimates the diameter of the spectrahedron. It samples points uniformly with coordinate directions
    /// hit and run, and returns the maximum distance between them.
    /// \tparam Point
    /// \param[in] numPoints The number of points to sample for the estimation
    /// \return An estimation of the diameter of the spectrahedron
    template<class Point>
    NT estimateDiameter(int const numPoints, Point const & interiorPoint) {
        typedef boost::mt19937 RNGType;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);

        std::list<Point> randPoints;

        // initialize random numbers generators
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(1, d);

        PrecomputedValues precomputedValues;
        VT p = interiorPoint.getCoefficients();

        // sample points with walk length set to 1
        for (int samplingNo=0 ; samplingNo<numPoints ; ++samplingNo) {
            // uniformly select a line parallel to an axis,
            // i.e. an indicator i s.t. x_i = 1
            int coordinate = uidist(rng);

            // get the distances we can travel from p
            // on the line p + t* e_coordinate
            // before reaching the boundary
            std::pair<NT, NT> distances = this->coordinateIntersection(p, coordinate, precomputedValues);

            // uniformly set the new point on the segment
            // defined by the intersection points
            NT lambda = urdist(rng);
            NT diff = distances.first + lambda * (distances.second - distances.first);

            p(coordinate - 1) = p(coordinate - 1) + diff;

            // update the precomputedValues, so we can skip
            // computations in the next call
            precomputedValues.computed_A = true;
            precomputedValues.A += diff * *(this->getLMI().getMatrix(coordinate));
            randPoints.push_back(Point(p));
        }

        // find maximum distance among points;
        NT maxDistance = 0;
        typename std::list<Point>::iterator itInner, itOuter = randPoints.begin();

        for (; itOuter!=randPoints.end() ; ++itOuter)
            for (itInner=itOuter ; itInner!=randPoints.end() ; ++itInner) {
                NT current = itOuter->distance(*itInner);
                if (current > maxDistance)
                    maxDistance = current;
            }

        return maxDistance;
    }

    /// Find out is lmi(current position) = mat is in the exterior of the spectrahedron
    /// \param mat a matrix where mat = lmi(current position)
    /// \return true if position is outside the spectrahedron
    bool isExterior(MT const & mat) {
        return !lmi.isNegativeSemidefinite(mat);
    }

    /// Find out is pos is in the exterior of the spectrahedron
    /// \param pos a vector
    /// \return true if pos is outside the spectrahedron
    bool isExterior(VT const & pos) {
        return !lmi.isNegativeSemidefinite(pos);
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
