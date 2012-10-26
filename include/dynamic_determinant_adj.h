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

#ifndef DYNAMIC_DETERMINANT_ADJ_H
#define DYNAMIC_DETERMINANT_ADJ_H

#include <ostream>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/assertions.h>
#include "indexed_point.h"
#include "sort_swap_adj.h"
#include "eigen_functions.h"
#include <vector>
#define tvector std::vector

namespace HeaDaCHe{

// This object stores the necessary data to compute a dynamic determinant.
// It can be constructed from the determinant and the inverse of a matrix A
// and the indices of the points that conform A, or by a range of d+1
// indexed points of dimension d.
template <class _NT, class _P>
class DeterminantData{
        public:
        typedef _NT                                     NT;
        //typedef _P																			Point;
        typedef tvector<NT>                             row_type;
        typedef tvector<row_type>                       matrix_type;
        private:
        typedef CGAL::Linear_algebraCd<NT>              LA;
        typedef typename LA::Matrix                     LAMatrix;
        typedef typename LA::Vector                     LAVector;

        public:
        // Default constructor. TODO: remove!
        DeterminantData(){};
        // Constructor from inverse and determinant of a matrix.
        DeterminantData(const matrix_type &m,
                        const NT &d,
                        const tvector<size_t> &c):
                _adjA(m),_detA(d),_colsA(c){};
        // Constructor that takes d+1 points in dimension d.
        template <class InputIterator>
        DeterminantData(InputIterator begin,InputIterator end){
                initialize(begin,end);
        };
				       
        // Initializes with d+1 points in dimension d.
				template <class InputIterator>
				void initialize(InputIterator begin,InputIterator end){
                size_t dimp1=std::distance(begin,end); // dimension+1
                LAMatrix m(dimp1);
                size_t col=0;
                _colsA.reserve(dimp1);                
                for(InputIterator it=begin;it!=end;++it){
                        for(size_t row=0;row<dimp1-1;++row){
                                m(col,row)=(*it)[row];
                                //std::cout<<m(col,row)<<" ";
                        }
                        m(col,dimp1-1)=1;
                        _colsA.push_back(it->index());
                        ++col;
                }
                //std::cout<<std::endl;
                LAMatrix adj = eigen_naive_adjoint(m);
                // Copy the adjoint matrix to our representation.
                tvector<NT> tmprow;
                tmprow.reserve(dimp1);
                _adjA.reserve(dimp1);
                for(size_t row=0;row<dimp1;++row){
                        tmprow.clear();
                        for(size_t col=0;col<dimp1;++col)	{												
                                tmprow.push_back(adj(col,row));
																//std::cout<<adj(col,row)<<" ";
                        }
                        _adjA.push_back(tmprow);
                }
                _detA=eigen_determinant(m);
                //std::cout<<_detA<<std::endl;
                
                //if the dimension is odd change the sign
                //to follow triangulation code 
                //TODO: check if needed
                //_detA = (dimp1-1)%2==1 ? _detA*NT(-1) : _detA;
                //std::cout <<"#det="<<_detA<<"\n";
        };
        const matrix_type& get_adjmatrix()const{return _adjA;};
        const NT& get_determinant()const{return _detA;};
        const tvector<size_t>& get_columns()const{return _colsA;};
        std::ostream& print_inverse(std::ostream &o)const{
                // The matrix is printed in rows, exactly as it is stored.
                o<<"[ ";
                for(size_t row=0;row<_adjA.size();++row){
                        o<<"[ ";
                        for(size_t col=0;col<_adjA[row].size();++col){
                                o<<_adjA[row][col]<<' ';
                        }
                        o<<"] ";
                }
                return o<<']';
        };
        std::ostream& print_columns (std::ostream &o)const{
                o<<"{ "<<std::flush;
                for(tvector<size_t>::iterator i=_colsA.begin();
                    i!=_colsA.end();
                    ++i)
                        o<<(*i)<<' '<<std::flush;
                return o<<'}'<<std::flush;
        };

				std::ostream& print_columns (std::ostream &o){
                o<<"{ "<<std::flush;
                for(tvector<size_t>::iterator i=_colsA.begin();
                    i!=_colsA.end();
                    ++i)
                        o<<(*i)<<' '<<std::flush;
                return o<<'}'<<std::flush;
        };

        private:
        matrix_type _adjA; // The adjoint of the matrix A.
        NT _detA; // The determinant of A.
        tvector<size_t> _colsA; // The indices of the points that
                                     // are the columns of A.
        //tvector<_P*> _pointsA; // pointers to the points of A
};


 
//using Sherman-Morisson, and Adjoint matrix
template <class _NT,class _P>
DeterminantData<_NT,_P>
ddeterminant(const DeterminantData<_NT,_P> &DD,
             const _P &Pold,
             const _P &Pnew){
        typedef _NT                                     NT;
        typedef DeterminantData<NT,_P>                  DetData;
        typedef typename DetData::matrix_type           matrix_type;
        const matrix_type &Aadj=DD.get_adjmatrix();
        const NT &detA=DD.get_determinant();
        CGAL_assertion_msg(detA!=0,"the stored determinant must not be zero");
        tvector<size_t> colsAp=DD.get_columns();
        size_t s=colsAp.size();//Pold.dimension()+1;
        size_t pos=0; // The column of A that will be replaced.
        while(colsAp[pos]!=Pold.index()){
                ++pos;
                CGAL_assertion(pos<s);
        }
        CGAL_assertion(colsAp[pos]==Pold.index());
        colsAp[pos]=Pnew.index();
        //
        // dif = u-(A)i
        tvector<NT> dif;
        dif.reserve(s);
        //std::cout<<Pnew.index()<<" "<<Pold.index()<<std::endl;
        for(size_t i=0;i<s-1;++i){
          dif.push_back(Pnew[i]-Pold[i]);
          //std::cout<<Pnew[i]<<" "<<Pold[i]<<std::endl;
        }       
        dif.push_back(0);
        // 				
        // A1 = Aadj*dif
        tvector<NT> A1;
        A1.reserve(s);
        for(size_t i=0;i<s;++i){
					NT sum=0;
					//std::cout<<i<<" "<<Aadj.size()<<std::endl;
					for(size_t j=0;j<s;++j){
						//std::cout<<i<<","<<j<<" "<<Aadj[i][j]<<std::endl;
						sum+=Aadj[i][j]*dif[j];
					}
					A1.push_back(sum);
        }
        //
        // detAp = detA + e_i^T * Aadj * dif
        NT detAp = detA + A1[pos];
				// 
				// Apadj means "A prime adjoint"
				// Apadj = A^{-1} - A1 * A^{-1}[pos]
				matrix_type Apadj; 
        tvector<NT> r;
        Apadj.reserve(s);
        r.reserve(s);
        NT div_exact_result;
        for(size_t i=0;i<s;++i){
					r.clear();
					for(size_t j=0;j<s;++j){
						/*
						mpz_divexact(div_exact_result.get_mpz_t(),
						        mpz_class(Aadj[i][j]*detAp-A1[i]*Aadj[pos][j]).get_mpz_t(),
						            detA.get_mpz_t());
						*/
						/*
						mpz_divexact(div_exact_result.mpz(),
						            (Aadj[i][j]*detAp-A1[i]*Aadj[pos][j]).mpz(),
						            detA.mpz());
									
						std::cout<<"\ndet from ADJ = ("<<Aadj[i][j]<<"*"<<detAp.mpz()<<"-"<<A1[i]<<"*"<< Aadj[pos][j]<<")/"<<detA<<" = "<<
						  (Aadj[i][j]*detAp-A1[i]*Aadj[pos][j])/detA<<"="
						  <<div_exact_result<<std::endl;
						*/
						/*
						CGAL_assertion(mpz_divisible_p(
						          (Aadj[i][j]*detAp-A1[i]*Aadj[pos][j]).mpz(),
						          detA.mpz()));
						/**/
						div_exact_result=(Aadj[i][j]*detAp-A1[i]*Aadj[pos][j])/detA;
						r.push_back(div_exact_result);
					}
					Apadj.push_back(r);
        }
        //std::cout << "s"<<s<<" indet"<<detAp<<std::endl;
        //detAp = (s-1)%2==1 ? detAp*NT(-1) : detAp;
        return DetData(Apadj,detAp,colsAp);
};

} // namespace HeaDaCHe

#endif // DYNAMIC_DETERMINANT_ADJ_H
// vim: ts=2:expandtab
