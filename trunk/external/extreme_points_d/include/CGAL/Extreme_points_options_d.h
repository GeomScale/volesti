// Copyright (c) 2010 ETH Zurich (Switzerland)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// 
// $URL: $
// $Id:  $
// 
//
// Author(s)     : Christian Helbling

#ifndef CGAL_EXTREME_POINTS_OPTIONS_D_H
#define CGAL_EXTREME_POINTS_OPTIONS_D_H

#include <CGAL/QP_models.h>

namespace CGAL {

enum Extreme_point_algorithm_d {
    EP_CHOOSE_APPROPRIATE,
    EP_SIMPLE,
    EP_DULA_HELGASON
};

class Extreme_points_options_d {
private:
    Extreme_point_algorithm_d algo_;
    Quadratic_program_options qp_options_;
    
public:
    // default constructor
    Extreme_points_options_d(Extreme_point_algorithm_d algo = EP_CHOOSE_APPROPRIATE,
                Quadratic_program_options qp_options = Quadratic_program_options())
      : algo_(algo), qp_options_(qp_options) {}
    
    // set/get algorithm
    // ------------------------
    Extreme_point_algorithm_d get_algorithm() const
    {
        return algo_;
    }
    
    void set_algorithm (Extreme_point_algorithm_d algo)
    {
        algo_ = algo;
    }
    
    // set/get qp_options
    // ------------------------
    Quadratic_program_options get_qp_options() const
    {
        return qp_options_;
    }
    
    void set_qp_options (Quadratic_program_options qp_options)
    {
        qp_options_ = qp_options;
    }
};


} //namespace CGAL

#endif // CGAL_EXTREME_POINTS_OPTIONS_D_H
