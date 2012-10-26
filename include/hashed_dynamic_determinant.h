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

#ifndef HASHED_DYNAMIC_DETERMINANT_H
#define HASHED_DYNAMIC_DETERMINANT_H

#include <vector>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
//#include "dynamic_determinant.h"
//#include "sort_swap.h"
#include "dynamic_determinant_adj.h"
#include "sort_swap_adj.h"
#include "eigen_functions.h"
#include <algorithm>


using namespace boost::tuples;

namespace HeaDaCHe{
	
template <class P>
struct compare_o:
public std::binary_function<P,P,bool>
 {
  bool operator() (const P &p1, const P &p2)const 
   { return p1.get<1>() < p2.get<1>();}
	};

template <class _R>
class HashedDynamicDeterminant{
	
	typedef _R															        R;
	typedef typename R::RT                          NT;
	typedef typename R::Point_d  										Point;
  typedef HeaDaCHe::DeterminantData<NT,Point>		  DD;
  typedef std::vector<size_t>                     Index;
  typedef boost::tuple<Point,size_t,DD>						point_DD;
  //typedef std::tuple<Point,size_t,DD>						point_DD;
	typedef boost::unordered_map<Index,point_DD>    DDMap;
	typedef typename DDMap::iterator								DDMap_iterator;
	typedef Index::const_iterator      							Index_it;
	
	typedef boost::tuple<Point,size_t,size_t> Point_idx; //a point with many indices
	
	public:
	
	// constructor for dynamic det table
	HashedDynamicDeterminant() {};
	
	// initialize DDmap
	template <class InputIterator>
	NT initialize(InputIterator first,InputIterator last){
		DD initial_dd(first,last);
		//std::cout<<"init_det="<<initial_dd.get_determinant()<<" ";
		update_det(first,last,initial_dd);
		return 0;
	}
	
	// compute the dynamic determinant
	template <class InputIterator>
	NT dynamic_determinant(InputIterator first,InputIterator last){
		// find the max indexed point (i.e. the newest) and its place
		//also create query idx
		Index query_idx;
		size_t new_point_idx=first->index();
		//the new point is the maximum indexed point
		Point new_point=*first;
		//std::cout<<"(";
		size_t new_point_loc_idx=0;
		for(InputIterator s2=first; s2!=last; ++s2 ){
				//std::cout<<s2->index()<<" ";
				query_idx.push_back(s2->index());
				if (s2->index()>new_point_idx){
					new_point_idx=s2->index();
					new_point_loc_idx=s2-first;
					new_point=*s2;
				}
		}
		//std::cout<<"){"<<new_point_loc_idx<<"}";
		//std::cout<<"max="<<new_point_idx<<" "<<std::flush;
		// construct the indices of all points except max
		// TODO: do it in the previous for loop
		Index subidx;
		//std::vector<FT> r;
		for(InputIterator it=first;it!=last;++it){
				if(*it!=new_point){        
						subidx.push_back(it->index());
						//r.push_back((*it)[it->dimension()-1]);
				}
		}
		//size_t max=new_point.index();
		//std::cout<<"max="<<max<<" "<<std::flush;
		//std::cout<<"subidx="<<subidx[0]<<","<<subidx[1]<<"\t"<<std::flush;
		//! sort the index and search
		Index idx_sorted = sort_swap(subidx).first;
		DDMap_iterator ddm_it =ddm.find(idx_sorted);
		//std::cout<<idx_sorted<<std::endl;
		
		//point_DD pDD=ddm[subidx];
		NT newdet(99999);
		if (ddm_it != ddm.end()){
			//std::cout<<"{DET_PRECOMPUTED}"<<std::flush;
			point_DD pDD=(*ddm_it).second;
		  //std::cout<<"det="<<pDD.second.get_determinant()<<" "<<std::flush;
			//std::cout<<"old_point:"<<pDD.first.index()<<" "<<std::flush;
			Index new_idx = pDD.get<2>().get_columns();
			new_idx[pDD.get<1>()]=new_point.index();
			
			//! COMPUTE DET
			DD new_dd=HeaDaCHe::ddeterminant(pDD.get<2>(),pDD.get<0>(),new_point);
			//std::cout<<"sort_dist="<<sort_distance(query_idx,new_idx)<<" ";	
			
			newdet = new_dd.get_determinant();
			newdet = sort_distance(query_idx,new_idx) ? (-1)*newdet : newdet;
			
			size_t dist=new_point_loc_idx>pDD.get<1>()?
			            new_point_loc_idx-pDD.get<1>():
			            pDD.get<1>()-new_point_loc_idx;
			//std::cout<<pDD.get<1>()<<" d={"<<dist<<"}";
			
			//! if the new det is not zero UPDATE
			if(newdet!=NT(0)){
				update_det(first,last,new_dd);
			}
		}else{ //compute it from scratch
			
			typedef typename R::LA          LA;
			int d = static_cast<int>(std::distance(first,last)) - 1;
			typename LA::Matrix M(d);
			InputIterator s = first;
      ++s;
      for( int j = 0; j < d; ++s, ++j )
          for( int i = 0; i < d; ++i )
              M(i,j) = s->cartesian(i) - first->cartesian(i);
      /*typedef typename R::LA          LA;
			int d = static_cast<int>(std::distance(first,last)) -1;
			typename LA::Matrix M(d+1);
			InputIterator s = first;
      //++s;
      for( int j = 0; j <= d; ++s, ++j ){
          for( int i = 0; i < d; ++i ){
              //M(i,j) = s->cartesian(i) - first->cartesian(i);
              M(i,j) = s->cartesian(i);
							std::cout<<M(i,j)<<" "<<std::flush;
          }
          M(d,j) = NT(1);
      }*/
      //std::cout<<"scratchscratchscratch"<<std::endl;
      newdet = eigen_determinant(M);
			//std::cout << "not found" << std::endl;
			/*
			DD initial_dd(first,last);
			newdet = initial_dd.get_determinant();
		  update_det(first,last,initial_dd);
		   */
		}
		return newdet;
	}
	
	// update DDMap 
	template <class InputIterator>
	void update_det(InputIterator first,
									InputIterator last,
									DD new_dd){		
		int d=std::distance(first,last);
		//d=D+1(D=dimension) consider all the possible combinations of D points
		//Index idx;
		std::vector<Point_idx> idxP;
		//for (InputIterator it=first; it!=last; ++it){
		//	size_t local_idx=std::distance(first,it);
		//	idxP.push_back(Point_idx(*it,it->index(),local_idx));
		//}
		//Index idx=new_dd.get_columns();
		size_t k=0;
		for (Index::const_iterator iit=(new_dd.get_columns()).begin(); 
		               iit!=(new_dd.get_columns()).end(); ++iit){
				for (InputIterator it=first; it!=last; ++it){
						if(it->index()==*iit){
								//size_t local_idx=std::distance(first,it);
								idxP.push_back(Point_idx(*it,it->index(),k++));
						}
				}
		}
		/*	
		for (typename std::vector<Point_idx>::const_iterator idx_it=idxP.begin(); 
		     idx_it!=idxP.end(); ++idx_it)
		     std::cout << idx_it->get<1>() << "(" 
										<< idx_it->get<2>()
										<<") "<<std::flush;
		
		*/
		std::sort(idxP.begin(),idxP.end(),compare_o<Point_idx>());
		/*
		for (typename std::vector<Point_idx>::const_iterator idx_it=idxP.begin(); 
		     idx_it!=idxP.end(); ++idx_it)
		     std::cout << idx_it->get<1>() << "(" 
										<< idx_it->get<2>()
										<<") "<<std::flush;
		*/
		size_t i=0;
		for (typename std::vector<Point_idx>::const_iterator idx_it=idxP.begin(); 
		     idx_it!=idxP.end(); ++idx_it){
			//std::cout << idx_it->get<1>() << "|"<<std::flush;
			//std::cout<<"["<<std::flush;
			Index subidx;
			
			for (size_t j=0; j<d; ++j){
				if (j!=i){
					//std::cout<<idxP[j].get<1>()<<","<<std::flush;
					subidx.push_back(idxP[j].get<1>());
				}
			}
			++i;
			//std::cout<<"]"<<std::flush;
			//if idx is not already in ddm insert it
			//TODO:don't check for the one you know it is in ddm
			if(ddm.find(subidx) == ddm.end()){
				//point_DD pDD((*(first+i)),i,new_dd);
				point_DD pDD(idx_it->get<0>(),idx_it->get<2>(),new_dd);
				ddm[subidx]=pDD;
			}
		}
		/*
		//Index sorted_idx = sort_swap(idx).first;
		//std::cout<<idx<<"-->"<<sorted_idx<<std::endl;
		for (int i=0; i<d; ++i){	
			Index subidx;
			std::cout<<"[";
			for (int j=0; j<d; ++j){
				if (j!=i){
					std::cout<<idx[j]<<" ";
					subidx.push_back(idx[j]);
				}
			}
			std::cout<<"]"<<std::flush;
			//if idx is not already in ddm insert it
			//TODO:don't check for the one you know it is in ddm
			if(ddm.find(idx) == ddm.end()){
				point_DD pDD((*(first+i)),i,new_dd);
				ddm[idx]=pDD;
			}
		}*/
	}
	
	// compute the dynamic determinant WITHOUT UPDATE
	template <class InputIterator>
	NT dynamic_determinant_only(InputIterator first,InputIterator last){
		// find the max indexed point (i.e. the newest) and its place
		//also create query idx
		Index query_idx;
		size_t new_point_idx=first->index();
		//the new point is the maximum indexed point
		Point new_point=*first;
		//std::cout<<"(";
		size_t new_point_loc_idx=0;
		for(InputIterator s2=first; s2!=last; ++s2 ){
				//std::cout<<s2->index()<<" ";
				query_idx.push_back(s2->index());
				if (s2->index()>new_point_idx){
					new_point_idx=s2->index();
					new_point_loc_idx=s2-first;
					new_point=*s2;
				}
		}
		//std::cout<<"){"<<new_point_loc_idx<<"}";
		//std::cout<<"max="<<new_point_idx<<" "<<std::flush;
		// construct the indices of all points except max
		// TODO: do it in the previous for loop
		Index subidx;
		//std::vector<FT> r;
		for(InputIterator it=first;it!=last;++it){
				if(*it!=new_point){        
						subidx.push_back(it->index());
						//r.push_back((*it)[it->dimension()-1]);
				}
		}
		//size_t max=new_point.index();
		//std::cout<<"max="<<max<<" "<<std::flush;
		//std::cout<<"subidx="<<subidx[0]<<","<<subidx[1]<<"\t"<<std::flush;
		//! sort the index and search
		Index idx_sorted = sort_swap(subidx).first;
		DDMap_iterator ddm_it =ddm.find(idx_sorted);
		//std::cout<<idx_sorted<<std::endl;
		
		//point_DD pDD=ddm[subidx];
		NT newdet(99999);
		if (ddm_it != ddm.end()){
			//std::cout<<"{DET_PRECOMPUTED}"<<std::flush;
			point_DD pDD=(*ddm_it).second;
		  //std::cout<<"det="<<pDD.second.get_determinant()<<" "<<std::flush;
			//std::cout<<"old_point:"<<pDD.first.index()<<" "<<std::flush;
			Index new_idx = pDD.get<2>().get_columns();
			new_idx[pDD.get<1>()]=new_point.index();
			
			//! COMPUTE DET
			DD new_dd=HeaDaCHe::ddeterminant(pDD.get<2>(),pDD.get<0>(),new_point);
			//std::cout<<"sort_dist="<<sort_distance(query_idx,new_idx)<<" ";	
			
			newdet = new_dd.get_determinant();
			newdet = sort_distance(query_idx,new_idx) ? (-1)*newdet : newdet;
			
			size_t dist=new_point_loc_idx>pDD.get<1>()?
			            new_point_loc_idx-pDD.get<1>():
			            pDD.get<1>()-new_point_loc_idx;
			//std::cout<<pDD.get<1>()<<" d={"<<dist<<"}";
			//std::cout<<newdet<<" ";
			return newdet;
			
		}else{ //compute it from scratch
			
			typedef typename R::LA          LA;
			int d = static_cast<int>(std::distance(first,last)) - 1;
			typename LA::Matrix M(d);
			InputIterator s = first;
      ++s;
      for( int j = 0; j < d; ++s, ++j )
          for( int i = 0; i < d; ++i )
              M(i,j) = s->cartesian(i) - first->cartesian(i);
      /*typedef typename R::LA          LA;
			int d = static_cast<int>(std::distance(first,last)) -1;
			typename LA::Matrix M(d+1);
			InputIterator s = first;
      //++s;
      for( int j = 0; j <= d; ++s, ++j ){
          for( int i = 0; i < d; ++i ){
              //M(i,j) = s->cartesian(i) - first->cartesian(i);
              M(i,j) = s->cartesian(i);
							std::cout<<M(i,j)<<" "<<std::flush;
          }
          M(d,j) = NT(1);
      }*/
      //std::cout<<"scratchscratchscratch"<<std::endl;
      newdet = eigen_determinant(M);
			//std::cout << "not found" << std::endl;
			/*
			DD initial_dd(first,last);
			newdet = initial_dd.get_determinant();
		  update_det(first,last,initial_dd);
		   */
		}
		return newdet;
	}
	
	// print DDMap
	void print(){
		std::cout<<"Dynamic Determinant Map:"<<std::endl;
		for(DDMap_iterator ddm_it=ddm.begin(); ddm_it!=ddm.end(); ++ddm_it){
			std::cout<<"["<<(*ddm_it).first<<"], "
							 <<(*ddm_it).second.get<0>().index()<<"("
							 <<(*ddm_it).second.get<1>()<<") --> "
			         <<(*ddm_it).second.get<2>().get_determinant()<<" ";
			(*ddm_it).second.get<2>().print_columns(std::cout);
			std::cout<<std::endl;
		}
		/*DDMap_iterator ddm_it=ddm.begin();
		std::vector<size_t> test;
		test.push_back(3);
		test.push_back(0);
		test.push_back(5);

		std::cout<<sort_distance((*ddm_it).second.get<2>().get_columns(),
									test)<<std::endl;
		*/
	}
	
        // get DDMap's size
        size_t get_size(){
                return ddm.size();
        }

	private:
	
	DDMap ddm;	
};

}//namespace HeaDaCHe

#endif // HASHED_DYNAMIC_DETERMINANT_H
