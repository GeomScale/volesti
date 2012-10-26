// Copyright (c) 2009 INRIA Sophia-Antipolis (France),
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
// $URL$
// $Id$
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_REGULAR_TRIANGULATION_H
#define CGAL_REGULAR_TRIANGULATION_H

#include <CGAL/Triangulation.h>
#include <CGAL/Dimension.h>
#include <CGAL/Default.h>

namespace CGAL {

template< typename RTTraits, typename TDS_ = Default >
class Regular_triangulation
: public Triangulation<RTTraits,
         typename Default::Get<TDS_, Triangulation_data_structure<
                             typename Ambient_dimension<typename RTTraits::Point_d>::type,
                             Triangulation_vertex<RTTraits>,
                             Triangulation_full_cell<RTTraits> >
                    >::type >
{
    typedef typename Ambient_dimension<typename RTTraits::Point_d>::type
                                                    Ambient_dimension_;
    typedef typename Default::Get<TDS_, Triangulation_data_structure<
                         Ambient_dimension_,
                         Triangulation_vertex<RTTraits>,
                         Triangulation_full_cell<RTTraits> >
                >::type                         TDS;
    typedef Triangulation<RTTraits, TDS>        Base;
    typedef Regular_triangulation<RTTraits, TDS_>    Self;

public:
    typedef Ambient_dimension_                  Ambient_dimension;
};

} //namespace CGAL

#endif CGAL_REGULAR_TRIANGULATION_H
