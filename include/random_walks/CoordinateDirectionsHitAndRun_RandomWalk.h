//
// Created by panagiotis on 2/28/20.
//

#ifndef VOLESTI_COORDINATEDIRECTIONSHITANDRUN_RANDOMWALK_H
#define VOLESTI_COORDINATEDIRECTIONSHITANDRUN_RANDOMWALK_H

/// The coordinate directions hit and run random walk, to sample a spectrahedron from the uniform distribution.
/// \tparam Point Point type to communicate with the rest library
/// \tparam MT Matrix type for internal computations
/// \tparam VT Vector type for internal computations
/// \tparam RNGType
template<class Point, typename MT, typename VT, typename RNGType>
class CoordinateDirectionsHitAndRun_RandomWalk {
public:

    /// The numeric type
    typedef typename Point::FT NT;
    /// The type of the spectrahedron
    typedef Spectrahedron<NT, MT, VT> SPECTRAHEDRON;
    /// Type for internal structure of class Spectrahedron
    typedef typename SPECTRAHEDRON::PrecomputedValues PrecomputedValues;

    /// The number of points to "burn", before keeping the following as a sample
    int walkLength;
    /// For generating random numbers
    RNGType rng;

    /// Constructor
    /// \param[in] walkLength The number of points to "burn", before keeping the following as a sample
    /// \param[in] rng For generating random numbers
    CoordinateDirectionsHitAndRun_RandomWalk(int const walkLength, RNGType const & rng ) : walkLength(walkLength), rng(rng) {}

    /// Samples uniformly random points in the spectrahedron
    /// \param[in] spectrahedron The spectrahedron to be sampled
    /// \param[in] interiorPoint A point in the interior of the spectrahedron
    /// \param[in] numPoints The number of points to sample
    /// \param[out] randPoints The list which will hold the samples
    void sample(SPECTRAHEDRON &spectrahedron, Point const & interiorPoint, const int numPoints, std::list<Point> & randPoints) {

        // store intermediate results between successive calls of methods
        // of the class spectrahedron, to avoid repeating computations
        PrecomputedValues precomputedValues;
        VT p = interiorPoint.getCoefficients();

        // sample #pointsNum points
        for (unsigned int i = 1; i <= numPoints; ++i) {
            // burn #walk_length points to get one sample
            for (unsigned int j = 0; j < walkLength; ++j) {
                getNextPoint(spectrahedron, p, precomputedValues);
            }

            // add the sample in the return list
            randPoints.push_back(Point(p));
        }

        // the data in preComputedValues may be out of date in the next call
        precomputedValues.resetFlags();
    }


    void getNextPoint(SPECTRAHEDRON & spectrahedron, VT &p, PrecomputedValues & precomputedValues) {
        int n = spectrahedron.dimension();

        // initialize random numbers generators
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(1, n);

        // uniformly select a line parallel to an axis,
        // i.e. an indicator i s.t. x_i = 1
        int coordinate = uidist(rng);

        // get the distances we can travel from p
        // on the line p + t* e_coordinate
        // before reaching the boundary
        std::pair<NT, NT> distances = spectrahedron.coordinateIntersection(p, coordinate, precomputedValues);

        // uniformly set the new point on the segment
        // defined by the intersection points
        NT lambda = urdist(rng);
        NT diff = distances.first + lambda*(distances.second - distances.first);

        p(coordinate-1) = p(coordinate-1) + diff;

        // update the precomputedValues, so we can skip
        // computations in the next call
        precomputedValues.computed_A = true;
        precomputedValues.A += diff * *(spectrahedron.getLMI().getMatrix(coordinate));
    }
};
#endif //VOLESTI_COORDINATEDIRECTIONSHITANDRUN_RANDOMWALK_H
