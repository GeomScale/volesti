// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Modified by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "LMI.h"
#include "chrono"


/// Among successive calls of this class methods, we may need to pass data
/// from one call to the next, to avoid repeating computations, or to efficiently update values
/// Warning: this struct assists in many methods; perhaps for different methods use different instances
template <typename NT, typename MT, typename VT>
struct PrecomputationOfValues {

    /// These flags indicate whether the corresponding matrices are computed
    /// if yes, we can use them and not compute them fro scratch
    // TODO: avoid the use of flags
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

    void set_mat_size(int const& m)
    {
        A.setZero(m, m);
        B.setZero(m, m);
        C.setZero(m, m);

        X.setZero(2*m, 2*m);
        Y.setZero(2*m, 2*m);

        eigenvector.setZero(m);
    }
};


/// This class manipulates a spectrahedron, described by a Linear Matrix Inequality i.e. LMI
/// \tparam Point Point Type
template<typename Point>
class Spectrahedron {
public:

    /// The numeric/matrix/vector types we use
    typedef Point                                             PointType;
    typedef typename Point::FT                                NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>              VT;

    double maxDouble = std::numeric_limits<double>::max();

    /// The type of a pair of NT
    typedef std::pair<NT, NT> pairNT;

    typedef PrecomputationOfValues<NT, MT, VT> _PrecomputationOfValues;

    _PrecomputationOfValues precomputedValues;

    EigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;

    /// The dimension of the spectrahedron
    unsigned int d;
    VT grad;
    std::pair<PointType, NT> _inner_ball;

    /// The linear matrix inequality that describes the spectrahedron
    LMI<NT, MT, VT> lmi;

    Spectrahedron() {}

    /// Creates a spectrahedron
    /// \param[in] lmi The linear matrix inequality that describes the spectrahedron
    Spectrahedron(const LMI<NT, MT, VT>& lmi) : lmi(lmi) {
        d = lmi.dimension();
        //grad.setZero(d);
        precomputedValues.resetFlags();
        precomputedValues.set_mat_size(lmi.sizeOfMatrices());
    }

    void set_interior_point(PointType const& r)
    {
        _inner_ball.first = r;
    }

    std::pair<PointType, NT> ComputeInnerBall() {
        NT radius = maxDouble;

        for (unsigned int i = 0; i < dimension(); ++i) {

            std::pair<NT, NT> min_max = coordinateIntersection(_inner_ball.first.getCoefficients(), i+1);

            if (min_max.first < radius) radius = min_max.first;
            if (-min_max.second < radius) radius = -min_max.second;
        }

        radius = radius / std::sqrt(NT(dimension()));
        _inner_ball.second = radius;

        return std::pair<PointType, NT>(_inner_ball.first, radius);
    }

    std::pair<Point,NT> InnerBall() const
    {
        return _inner_ball;
    }

    /// Construct the quadratic eigenvalue problem \[At^2 + Bt + C \] for positive_intersect.
    /// A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// \param[in] a Input vector
    /// \param[in] b Input vector
    /// \param[in] c Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPositiveQuadIntersection(const VT& a, const VT& b, const VT& c) {

        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_A) {
            lmi.evaluateWithoutA0(a, precomputedValues.A, true);
        }

        if (!precomputedValues.computed_C) {
            lmi.evaluate(c, precomputedValues.C, true);
        }

        // compute Matrix B
        lmi.evaluateWithoutA0(b, precomputedValues.B, true);
    }

    /// Construct the generalized eigenvalue problem \[Bt + C \] for positive_intersect.
    /// \param[in] p Input vector
    /// \param[in] v Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPositiveLinearIntersection(const VT& p, const VT& v) {
        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_C) {
            lmi.evaluate(p, precomputedValues.C);
        }

        // compute Matrix B
        lmi.evaluateWithoutA0(v, precomputedValues.B);
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
    /// \returns The distance d
    NT positiveQuadIntersection(VT const & a, VT const & b, VT const & c) {
        unsigned int matrixDim = lmi.sizeOfMatrices();

        // create matrices A, B, C
        createMatricesForPositiveQuadIntersection(a, b, c);

        // get the minimum positive eigenvalue of At^2 + Bt + C
        NT distance = EigenvaluesProblem.minPosQuadraticEigenvalue(precomputedValues.A, precomputedValues.B,
                                                                   precomputedValues.C, precomputedValues.X,
                                                                   precomputedValues.Y,
                                                                   precomputedValues.eigenvector,
                                                                   precomputedValues.computed_XY);
        return distance;
    }


    NT positiveLinearIntersection(VT const & p, VT const & v)
    {
        createMatricesForPositiveLinearIntersection(p, v);
        NT distance = EigenvaluesProblem.minPosLinearEigenvalue(precomputedValues.C, precomputedValues.B,
                                                                precomputedValues.eigenvector);
        return distance;
    }

    /// Computes the distance d one must travel on the line a + tb,
    /// assuming we start at t=0 and that b has zero everywhere and 1 in its i-th coordinate.
    /// We must solve the generalized eigenvalue problem A+tB, where A = lmi(a) and B=(lmi) - A0 = A_i
    /// If the flag precomputedValues,computed_A is true, the matrix A is not computed.
    /// \param[in] a Input vector
    /// \param[in] coordinate Indicator of the i-th coordinate, 1 <= coordinate <= dimension
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT coordinateIntersection(VT const & a, int const coordinate) {

        // prepare the generalized eigenvalue problem A+lB
        // we may not have to compute A!
        if (!precomputedValues.computed_A)
            lmi.evaluate(a, precomputedValues.A);

        return EigenvaluesProblem.symGeneralizedProblem(precomputedValues.A, *(lmi.getMatrix(coordinate)));
    }

    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int const& rand_coord,
                                          VT&)
    {
        return coordinateIntersection(r.getCoefficients(), rand_coord);
    }

    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(PointType &r,
                                          PointType&,
                                          unsigned int const& rand_coord,
                                          unsigned int&,
                                          VT&)
    {
        return coordinateIntersection(r.getCoefficients(), rand_coord);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v)
    {
        NT pos_inter = positiveLinearIntersection(r.getCoefficients(), v.getCoefficients());
        return std::pair<NT, int> (pos_inter, -1);
    }

    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT& ,
                                               NT const&) {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT& ,
                                               NT const&,
                                               update_parameters&)
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT&,
                                               NT const&,
                                               MT const&,
                                               update_parameters& )
    {
        return line_positive_intersect(r, v);
    }

    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(PointType const& r,
                                                     PointType const& v,
                                                     VT&,
                                                     VT&,
                                                     update_parameters&)
    {
        return line_positive_intersect(r, v);
    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(PointType const& r,
                                               PointType const& v,
                                               VT&,
                                               VT&)
    {
        return line_positive_intersect(r, v);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(PointType const& r, PointType const& v)
    {
        NT pos_inter = positiveLinearIntersection(r.getCoefficients(), v.getCoefficients());
        NT neg_inter = -positiveLinearIntersection(r.getCoefficients(), NT(-1)*v.getCoefficients());

        return std::make_pair(pos_inter, neg_inter);
    }


    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT&,
                                    VT&)
    {
        return line_intersect(r, v);
    }

    std::pair<NT,NT> line_intersect(PointType const& r,
                                    PointType const& v,
                                    VT&,
                                    VT&,
                                    NT&)
    {
        return line_intersect(r, v);
    }

    void update_position_internal(NT &t){
        precomputedValues.C += t * precomputedValues.B;
        precomputedValues.computed_C = true;
    }

    MT get_mat() const
    {
        return MT::Identity(lmi.dimension(), lmi.dimension());
    }
    
    bool is_normalized ()
    {
        return true;
    }

    void normalize() {}

    void resetFlags()
    {
        precomputedValues.resetFlags();
    }

    void set_flags(bool bool_flag)
    {
        precomputedValues.computed_A = bool_flag;
        precomputedValues.computed_C = bool_flag;
        precomputedValues.computed_XY = bool_flag;
    }

    MT get_C() const
    {
        return precomputedValues.C;
    }

    void update_C(NT const& lambda)
    {
        precomputedValues.C += (lambda * lambda) * precomputedValues.A + lambda * precomputedValues.B;
    }

    // return the number of facets
    int num_of_hyperplanes() const
    {
        return 0;
    }

    void shift(VT e) {
        MT A0 = getLMI().get_A0();
        std::vector<MT> matrices = getLMI().getMatrices();

        int d = matrices.size();

        for (int i = 1; i < d; ++i) {
            A0 = A0 + e(i-1)*matrices[i];
        }

        lmi.set_A0(A0);

        _inner_ball.first = PointType(dimension());
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] r A point on the boundary of the spectrahedron
    /// \param[in] v The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    template <typename update_parameters>
    void compute_reflection(PointType &v, PointType const& r, update_parameters& ) const
    {
        VT grad(d);
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
    template <typename RNGType>
    NT estimateDiameter(int const numPoints, PointType const & interiorPoint, RNGType &rng) {

        std::list<Point> randPoints;

        precomputedValues.computed_A = false;
        VT p = interiorPoint.getCoefficients();

        // sample points with walk length set to 1
        for (int samplingNo=0 ; samplingNo<numPoints ; ++samplingNo) {
            // uniformly select a line parallel to an axis,
            // i.e. an indicator i s.t. x_i = 1
            int coordinate = rng.sample_uidist() + 1;

            // get the distances we can travel from p
            // on the line p + t* e_coordinate
            // before reaching the boundary
            std::pair<NT, NT> distances = this->coordinateIntersection(p, coordinate);

            // uniformly set the new point on the segment
            // defined by the intersection points
            NT lambda = rng.sample_urdist();
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

    int is_in(PointType const& p, NT tol=NT(0)) const
    {
        if (isExterior(p.getCoefficients())) {
            return 0;
        }
        return -1;
    }

    /// Find out is lmi(current position) = mat is in the exterior of the spectrahedron
    /// \param mat a matrix where mat = lmi(current position)
    /// \return true if position is outside the spectrahedron
    bool isExterior(MT const & mat) const {
        return !lmi.isNegativeSemidefinite(mat);
    }

    /// Find out is pos is in the exterior of the spectrahedron
    /// \param pos a vector
    /// \return true if pos is outside the spectrahedron
    bool isExterior(VT const & pos) const {
        return !lmi.isNegativeSemidefinite(pos);
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
