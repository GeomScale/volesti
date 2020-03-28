// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "LMI.h"
#include "../matrices/EigenvaluesProblems.h"

/// This class manipulates a spectrahedron, described by a LMI
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class Spectrahedron {
public:

    /// The type of a pair of NT
    typedef std::pair<NT, NT> pairNT;

    /// The dimension of the spectrahedron
    unsigned int d;

    /// The linear matrix inequality that describes the spectrahedron
    LMI<NT, MT, VT> lmi;

    /// Among successive calls of this class methods, we may need to pass data
    /// from one call to the next, to avoid repeating computations, or to efficiently update values
    /// Warning: this struct assists in many methods; perhaps for different methods use different instances
    struct PrecomputedValues {

        /// These flags indicate whether the corresponding matrices are computed
        /// if yes, we can use them and not compute them fro scratch
        bool computed_A;

        /// The matrices the method positiveIntersection receives from its previous call
        /// if the flag first_positive_intersection is true.
        /// Matrix A is also used in coordinateIntersection
        MT A;

        /// The eigenvector corresponding to the eigenvalue we kept
        /// from an eigenvalue problem we solved
        VT eigenvector;

        /// Sets all flags to false
        void resetFlags() {
            computed_A = false;
        }

        /// constructor - set defaule values
        PrecomputedValues() {
            computed_A = false;
        }
    };


    Spectrahedron() {}

    /// Creates a spectrahedron
    /// \param[in] lmi The linear matrix inequality that describes the spectrahedron
    Spectrahedron(const LMI<NT, MT, VT>& lmi) : lmi(lmi) {
        d = lmi.dimension();
    }

    /// \return The dimension of the spectrahedron
    unsigned int dimension() const {
        return d;
    }

    /// \return The LMI describing this spectrahedron
    LMI<NT, MT, VT> getLMI() const {
        return lmi;
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

    /// Computes the distance d one must travel on the line a + tb,
    /// assuming we start at t=0
    /// We must solve the generalized eigenvalue problem A+tB, where A = lmi(a) and B=(lmi) - A0
    /// If the flag precomputedValues,computed_A is true, the matrix A is not computed.
    /// \param[in] a Current position
    /// \param[in] b The direction we will follow
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT intersection(VT const & a, VT const & b, PrecomputedValues& precomputedValues) {

        // prepare the generalized eigenvalue problem A+lB
        // we may not have to compute A!
        if (!precomputedValues.computed_A)
            lmi.evaluate(a, precomputedValues.A);

        MT B;
        lmi.evaluateWithoutA0(b, B);

        EigenvaluesProblems<NT, MT, VT> eigenvaluesProblems;
        return eigenvaluesProblems.symGeneralizedProblem(precomputedValues.A, B);
    }

    /// Computes the distance d one must travel on the line a + tb,
    /// assuming we start at t=0 and we increase t till we reach the boundary
    /// We must solve the generalized eigenvalue problem A+tB, where A = lmi(a) and B=(lmi) - A0
    /// If the flag precomputedValues,computed_A is true, the matrix A is not computed.
    /// It also stores in precomputedValues the eigenvector from the eigenvalue problem we solved.
    /// It can be used to compute the reflection at the intersection point.
    /// \param[in] a Current position
    /// \param[in] b The direction we will follow
    /// \return The distance till we reach the boundary
    NT positiveIntersection(VT const & a, VT const & b, PrecomputedValues& precomputedValues) {

        // prepare the generalized eigenvalue problem A+lB
        // we may not have to compute A!
        if (!precomputedValues.computed_A)
            lmi.evaluate(a, precomputedValues.A);

        MT B;
        lmi.evaluateWithoutA0(b, B);

        EigenvaluesProblems<NT, MT, VT> eigenvaluesProblems;
        return eigenvaluesProblems.minPosSymGeneralizedProblem(precomputedValues.A, B, precomputedValues.eigenvector);
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
std::cout << grad;
        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * incomingDirection.dot(grad);
        reflectedDirection = incomingDirection - dot * grad;
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
