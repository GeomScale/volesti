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

#ifndef SORT_SWAP_ADJ_H
#define SORT_SWAP_ADJ_H

#include <cassert>
#include <boost/version.hpp>
#if BOOST_VERSION < 104300
#include <boost/detail/algorithm.hpp>
#else
#include <boost/range/algorithm_ext/is_sorted.hpp>
#endif

// This function takes as parameter an index vector v (usually, a
// std::vector<size_t>) and returns a pair, containing the sorted vector s
// and a boolean b. b is true iff s can be obtained from v with an even
// number of swaps.
// TODO: This function is naively implemented, it sorts using a bubble sort
// and computes the number of swaps on the fly. It will be better to
// implement a subquadratic algorithm that counts the swaps.
template <class Index>
std::pair<Index,bool> sort_swap(const Index &v){
        typedef typename Index::value_type                      elt_t;
        Index s(v); // The sorted vector.
        size_t swaps=0; // The number of swaps used.
        elt_t tmp; // The temporary for swaps.
        size_t min; // The temporary to store the minimum element.
        // Bubble sort, O(n^2) but easy to implement.
        for(size_t i=0;i+1<s.size();++i){
                min=i;
                // Find the minimum in v[i+1...].
                for(size_t j=i+1;j<s.size();++j)
                        if(s[j]<s[min])
                                min=j;
                // Swap s[i] and s[min] and update swaps if necessary.
                if(min!=i){
                        tmp=s[min];
                        s[min]=s[i];
                        s[i]=tmp;
                        ++swaps;
                }
        }
        assert(boost::is_sorted(s));
        assert(s.size()==v.size());
        return std::make_pair(s,!(swaps%2));
}

template <class Index>
bool sort_distance(const Index &v, const Index &w){
        // we compute the sum of the distances of v,w
        // from the sorted vector (^ is XOR)
        // TODO: do it more efficiently
        //e.g. distfromsorted(132)=1,distfromsorted(231)=2
        // dist(132,231)=1< distfromsorted(132)+distfromsorted(231)
        //std::cout<<"\n"<<"ind1="<<v<<" ";
        //std::cout<<"ind2="<<w<<" (";
        //std::cout<<"d1="<<!sort_swap(v).second<<" ";
        //std::cout<<"d2="<<!sort_swap(w).second<<")->{"
        //         << !(sort_swap(v).second ^ !sort_swap(w).second) <<"}\n";
        bool st = !sort_swap(v).second ^ !sort_swap(w).second;
        //if the dimension is odd change the sign
        //to follow triangulation code 
        int dim = v.size()-1;
        //std::cout<<"OLD:"<<st;
        st = dim%2==1 ? !st : st;
        //std::cout<<"NEW:"<<st; 
        return st;
}

#endif // SORT_SWAP_ADJ_H
// vim: ts=2
