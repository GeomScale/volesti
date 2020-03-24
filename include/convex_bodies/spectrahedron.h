//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "LMI.h"

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

};

#endif //VOLESTI_SPECTRAHEDRON_H
