// Copyright 2012 National and Kapodistrian University of Athens, Greece.
//
// This file is part of HeaDaCHe.
//
// HeaDaCHe is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// HeaDaCHe is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDaCHe,
// see <http://www.gnu.org/licenses/>.

#ifndef HEADACHE_H
#define HEADACHE_H

#include "dynamic_determinant.h"
#include "hashed_dynamic_determinant.h"
#include "indexed_point.h"
#include <CGAL/Triangulation_vertex.h>
#include <CGAL/Triangulation_full_cell.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Triangulation.h>

namespace HeaDaCHe{

template <class TrTraits>
class Triangulation:
public CGAL::Triangulation<
                TrTraits,
                CGAL::Triangulation_data_structure<
                        CGAL::Dynamic_dimension_tag,
                        CGAL::Triangulation_vertex<TrTraits>,
                        CGAL::Triangulation_full_cell<
                                TrTraits,
                                HeaDaCHe::DeterminantData
                                        <typename TrTraits::FT> > > >{

        typedef TrTraits                                        K;
        typedef CGAL::Triangulation_vertex<K>                   V;
        typedef HeaDaCHe::DeterminantData<typename K::FT>       DD;
        typedef CGAL::Triangulation_full_cell<K,DD>             FC;
        typedef CGAL::Triangulation_data_structure
                        <CGAL::Dynamic_dimension_tag,V,FC>      TDS;
        typedef CGAL::Triangulation<K,TDS>                      Base;

        public:
        Triangulation(const int dim,const K tt=K()):Base(dim,tt){};
        // TODO: implement copy constructor.
};

template <class _P>
struct point_container:public tvector<IndexedPoint<_P> >{
        typedef _P                                      Base;
        typedef IndexedPoint<Base>                      IPoint;

        template <class InputIterator>
        point_container(InputIterator begin,InputIterator end){
                size_t s=std::distance(begin,end);
                this->reserve(s);
                for(InputIterator i=begin;i!=end;++i)
                        this->push_back(IPoint(*i));
        }
};

} // namespace HeaDaCHe

#endif // HEADACHE_H
