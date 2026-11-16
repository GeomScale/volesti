#ifndef ROUNDING_HPP
#define ROUNDING_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <tuple>

#include "core_types.hpp"
#include "inscribed_ellipsoid_rounding.hpp"

template <int ellipsoid_type, typename MT, typename VT, typename NT,
          typename PolytopeType, typename PointType>
void run_rounding(PolytopeType& Polytope, PointType& center,
                  MT& T, VT& shift, NT& round_val)
{
    auto result = inscribed_ellipsoid_rounding<MT, VT, NT, PolytopeType, PointType,
                                               ellipsoid_type>(Polytope, center);
    T         = std::get<0>(result);
    shift     = std::get<1>(result);
    round_val = std::get<2>(result);
}

template <typename MT, typename VT, typename NT, typename PolytopeType, typename PointType>
void apply_polytope_rounding(const std::string& method,
                             PolytopeType& Polytope,
                             PointType& center,
                             MT& T,
                             VT& shift,
                             NT& round_val,
                             bool show_console_logs)
{
    if (show_console_logs) {
        std::cout << "[ROUNDING] Rounding is enabled. Applying " << method << " rounding...\n";
    }

    if (method == "max_ellipsoid")
        run_rounding<EllipsoidType::MAX_ELLIPSOID>(Polytope, center, T, shift, round_val);
    else if (method == "log_barrier")
        run_rounding<EllipsoidType::LOG_BARRIER>(Polytope, center, T, shift, round_val);
    else if (method == "vaidya_barrier")
        run_rounding<EllipsoidType::VAIDYA_BARRIER>(Polytope, center, T, shift, round_val);
    else if (method == "volumetric_barrier")
        run_rounding<EllipsoidType::VOLUMETRIC_BARRIER>(Polytope, center, T, shift, round_val);
    else
        throw std::runtime_error("Unknown rounding method: " + method);

    // Rounding shifts the polytope to the origin, so the new center is the zero vector.
    center = PointType(VT::Zero(Polytope.dimension()));
    if (show_console_logs) {
        std::cout << "[ROUNDING] Rounding complete. Round value: " << round_val << "\n\n";
    }
}

#endif // ROUNDING_HPP