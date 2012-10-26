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

#ifndef INDEXED_POINT_H
#define INDEXED_POINT_H

// An indexed point is a point that carries an index. An object of this
// type inherits from a point class and adds a data member, that is
// represented as a global variable (or thread-local, in multi-threaded
// contexts).

#include <CGAL/basic.h>

namespace HeaDaCHe {

#ifdef CGAL_HAS_THREADS
#include <boost/thread/tss.hpp>
static boost::thread_specific_ptr<std::size_t> max_index;
#else
static std::size_t max_index=0;
#endif

template <class _P>
class IndexedPoint:public _P{
        public:
        typedef _P                                      Base;

        public:
        IndexedPoint(const IndexedPoint &p):idx(p.index()),Base(p){std::cout<<1<<std::endl;}
        IndexedPoint(const Base &p):Base(p){
#ifdef CGAL_HAS_THREADS
                if(max_index.get()==NULL)
                        max_index.reset(new std::size_t(0));
                idx=*max_index;
                max_index.reset(new std::size_t(*max_index+1));
                std::cout << "CONSTRUCTOR: " << *max_index-1 << " created\n";
#else
								std::cout << "nothreadCONSTRUCTOR: " << max_index << " created\n";
                idx=max_index++;
#endif
        };
        std::size_t index()const{return idx;};
        
        private:
        std::size_t idx;
};

} // namespace HeaDaCHe

#endif // INDEXED_POINT_H
