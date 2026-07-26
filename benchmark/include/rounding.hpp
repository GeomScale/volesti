#ifndef ROUNDING_HPP
#define ROUNDING_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <tuple>

#include "core_types.hpp"
#include "inscribed_ellipsoid_rounding.hpp"

template <typename MT, typename VT, typename NT, typename PolytopeType, typename PointType>
void apply_polytope_rounding(const std::string& method, 
                             PolytopeType& Polytope, 
                             PointType& center, 
                             MT& T, 
                             VT& shift, 
                             NT& round_val) {
    
    std::cout << "[ROUNDING] Rounding is enabled. Applying " << method << " rounding...\n";

    if (method == "max_ellipsoid") {
        auto rounding_result = inscribed_ellipsoid_rounding<MT, VT, NT>(Polytope, center);
        T = std::get<0>(rounding_result);
        shift = std::get<1>(rounding_result);
        round_val = std::get<2>(rounding_result);
    } 
    else if (method == "log_barrier") {
        auto rounding_result = inscribed_ellipsoid_rounding<MT, VT, NT, decltype(Polytope), decltype(center), 2>(Polytope, center);
        T = std::get<0>(rounding_result);
        shift = std::get<1>(rounding_result);
        round_val = std::get<2>(rounding_result);
    } 
    else if (method == "vaidya_barrier") {
        auto rounding_result = inscribed_ellipsoid_rounding<MT, VT, NT, decltype(Polytope), decltype(center), 3>(Polytope, center);
        T = std::get<0>(rounding_result);
        shift = std::get<1>(rounding_result);
        round_val = std::get<2>(rounding_result);
    } 
    else if (method == "volumetric_barrier") {
        auto rounding_result = inscribed_ellipsoid_rounding<MT, VT, NT, decltype(Polytope), decltype(center), 4>(Polytope, center);
        T = std::get<0>(rounding_result);
        shift = std::get<1>(rounding_result);
        round_val = std::get<2>(rounding_result);
    } 
    else {
        throw std::runtime_error("Unknown rounding method: " + method);
    }
    
    // Since rounding shifts the polytope to the origin we will use a zero vector as the center.
    center = PointType(VT::Zero(Polytope.dimension()));
    
    std::cout << "[ROUNDING] Rounding complete. Round value: " << round_val << "\n\n";
}

#endif // ROUNDING_HELPER_HPP