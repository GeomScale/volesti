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

#ifndef CGAL_EXTREME_POINTS_TRAITS_D_H
#define CGAL_EXTREME_POINTS_TRAITS_D_H

#include <CGAL/basic_classes.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

#include <boost/functional.hpp>
#include <boost/function.hpp>
#include <boost/iterator.hpp>

namespace CGAL {

// this one is not really useful (except the Kernel just happens to fit in)
// but the specialized versions are..
template <class P>
class Extreme_points_traits_d {
    // example which would work for Point_d
public:
    typedef P                                           Point;
    typedef CGAL::Kernel_traits<Point>                  KTraits;
    typedef typename KTraits::Kernel                    Kernel;
    
    typedef typename Kernel::RT                         RT;
    typedef typename Kernel::Less_lexicographically_d   Less_lexicographically;
    struct Homogeneous_begin {
        typedef typename Point::Homogeneous_const_iterator result_type;
        result_type operator() (const Point &p) const {
            return p.homogeneous_begin();
        }
    };
};

template <class Kernel>
class Extreme_points_traits_d<Point_d<Kernel> >
{
public:  
    typedef Point_d<Kernel>                             Point;
    typedef typename Kernel::RT                         RT;
    typedef typename Kernel::Less_lexicographically_d   Less_lexicographically;
    struct Homogeneous_begin {
        typedef typename Point::Homogeneous_const_iterator result_type;
        result_type operator() (const Point &p) const {
            return p.homogeneous_begin();
        }
    };
};

template <class Kernel>
class Extreme_points_traits_d<Point_2<Kernel> > 
{
public:  
    typedef Point_2<Kernel>                             Point;
    typedef typename Kernel::RT                         RT;
    typedef typename Kernel::Less_xy_2                  Less_lexicographically;
    
    struct Homogeneous_begin {
        typedef boost::function<RT (int)> Function;
        typedef boost::counting_iterator<int> CI;
        typedef boost::transform_iterator<Function, CI> result_type;
        result_type operator() (const Point &p) const {
            return result_type(CI(0),
                               boost::bind1st(
                                   std::mem_fun(&Point::homogeneous), &p)
                              );
        }
    };
};

template <class Kernel>
class Extreme_points_traits_d<Point_3<Kernel> >
{
public:  
    typedef Point_3<Kernel>                             Point;
    typedef typename Kernel::RT                         RT;
    typedef typename Kernel::Less_xyz_3                 Less_lexicographically;
    
    struct Homogeneous_begin {
        typedef boost::function<RT (int)> Function;
        typedef boost::counting_iterator<int> CI;
        typedef boost::transform_iterator<Function, CI> result_type;
        result_type operator() (const Point &p) const {
            return result_type(CI(0),
                               boost::bind1st(
                                   std::mem_fun(&Point::homogeneous), &p)
                              );
        }
    };
};


} // namespace CGAL

#endif // CGAL_EXTREME_POINTS_TRAITS_D_H
