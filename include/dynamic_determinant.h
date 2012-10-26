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

#ifndef DYNAMIC_DETERMINANT_H
#define DYNAMIC_DETERMINANT_H

#include <ostream>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/assertions.h>
#include "indexed_point.h"
#include "sort_swap.h"
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
                _invA(m),_detA(d),_colsA(c){};
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
                        }
                        m(col,dimp1-1)=1;
                        _colsA.push_back(it->index());
                        ++col;
                }
                LAMatrix tmpi;
                LAVector tmpv;
                // The function LA::inverse performs Gaussian elimination
                // on the input matrix, and computes the determinant in the
                // process when the determinant is nonzero. If it is zero,
                // it doesn't compute the inverse, but a nontrivial
                // solution to the degenerate system (which we don't need).
                if(!LA::inverse(m,tmpi,_detA,tmpv)) // false iff _detA==0
                        return;
                //std::cout<<"$det="<<_detA<<" ";
                //_detA=LA::determinant(m);
                //std::cout<<"%det="<<_detA<<" ";
         
                // Copy the inverse matrix to our representation.
                tvector<NT> tmprow;
                tmprow.reserve(dimp1);
                _invA.reserve(dimp1);
                for(size_t row=0;row<dimp1;++row){
                        tmprow.clear();
                        for(size_t col=0;col<dimp1;++col)
                                tmprow.push_back(tmpi(col,row)/_detA);
                        _invA.push_back(tmprow);
                }
                _detA=eigen_determinant(m);
                //if the dimension is odd change the sign
                //to follow triangulation code 
                //TODO: check if needed
                _detA = (dimp1-1)%2==1 ? _detA*NT(-1) : _detA;
                //std::cout <<"#det="<<_detA<<"\n";
        };
        
        
        //void set_invmatrix(const matrix_type &m){_invA=m;};
        //void set_determinant(const NT &d){_detA=d;};
        const matrix_type& get_invmatrix()const{return _invA;};
        const NT& get_determinant()const{return _detA;};
        const tvector<size_t>& get_columns()const{return _colsA;};
        std::ostream& print_inverse(std::ostream &o)const{
                // The matrix is printed in rows, exactly as it is stored.
                o<<"[ ";
                for(size_t row=0;row<_invA.size();++row){
                        o<<"[ ";
                        for(size_t col=0;col<_invA[row].size();++col){
                                o<<_invA[row][col]<<' ';
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
        matrix_type _invA; // The inverse of the matrix A.
        NT _detA; // The determinant of A.
        tvector<size_t> _colsA; // The indices of the points that
                                     // are the columns of A.
        //tvector<_P*> _pointsA; // pointers to the points of A
};

// This function receives and object containing the dynamic determinant
// information for a matrix A and points Pold and Pnew. The point Pold will
// be replaced by the point Pnew in A, to form a new matrix A'. The
// function returns an object containing the dynamic determinant
// information for A'. _NT must be a field because this function uses
// divisions. _P is a point type, such as CGAL::Cartesian<_NT>::Point_d.
template <class _NT,class _P>
DeterminantData<_NT,_P>
ddeterminant_old(const DeterminantData<_NT,_P> &DD,
             const _P &Pold,
             const _P &Pnew){
        typedef _NT                                     NT;
        typedef DeterminantData<NT,_P>                  DetData;
        typedef typename DetData::matrix_type           matrix_type;
        // Compute first determinant(A').
        const matrix_type &Ai=DD.get_invmatrix();
        const NT &dA=DD.get_determinant();
        CGAL_assertion_msg(dA!=0,"the stored determinant must not be zero");
        tvector<size_t> colsAp=DD.get_columns();
        size_t s=colsAp.size();//Pold.dimension()+1;
        size_t pos=0; // The column of A that will be replaced.
        while(colsAp[pos]!=Pold.index()){
                ++pos;
                CGAL_assertion(pos<s);
        }
        CGAL_assertion(colsAp[pos]==Pold.index());
        colsAp[pos]=Pnew.index();
        tvector<NT> b,dif;
        b.reserve(s);
        dif.reserve(s);
        for(size_t i=0;i<s-1;++i)
                dif.push_back(Pnew[i]-Pold[i]);
        dif.push_back(0);
        for(size_t i=0;i<s;++i){
                NT temp=0;
                for(size_t j=0;j<s;++j)
                        temp+=Ai[i][j]*dif[j];
                b.push_back(temp);
        }
        NT bposp1=b[pos]+1; // This will be used a few times.

        // With this data, we can compute the determinant d=(1+b[pos])*dA.
        // If d==0, the inverse matrix cannot, and should not, be computed.
        // In that case, we return an empty matrix and d=0;
        matrix_type Api; // Api means "A prime inverse".
        if(bposp1==0)
                return DetData(Api,bposp1,colsAp);

        // Compute the inverse matrix (A')^(-1)=B^(-1)*A^(-1).
        Api.reserve(s);
        tvector<NT> r;
        r.reserve(s);
        for(size_t row=0;row<s;++row){ // row is the row of Api
                r.clear();
                if(row==pos){
                        for(size_t col=0;col<s;++col){ // column of Api
                                // row=pos is the easy case, Api[pos][col]=
                                r.push_back(Ai[pos][col]/bposp1);
                        }
                        Api.push_back(r);
                        ++row;
                        if(row==s)
                                break;
                        r.clear();
                }
                for(size_t col=0;col<s;++col){ // column of Api
                        // Api[row][col]=
                        r.push_back(Ai[row][col]-(b[row]*Ai[pos][col])/bposp1);
                }
                Api.push_back(r);
        }

        return DetData(Api,bposp1*dA,colsAp);
};

//using Sherman-Morisson with inverse matrix
template <class _NT,class _P>
DeterminantData<_NT,_P>
ddeterminant(const DeterminantData<_NT,_P> &DD,
             const _P &Pold,
             const _P &Pnew){
        typedef _NT                                     NT;
        typedef DeterminantData<NT,_P>                  DetData;
        typedef typename DetData::matrix_type           matrix_type;
        const matrix_type &Ai=DD.get_invmatrix();
        const NT &dA=DD.get_determinant();
        CGAL_assertion_msg(dA!=0,"the stored determinant must not be zero");
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
        for(size_t i=0;i<s-1;++i)
          dif.push_back(Pnew[i]-Pold[i]);       
        dif.push_back(0);
        // 				
        // A1 = A^{-1}*dif
        tvector<NT> A1;
        A1.reserve(s);
        for(size_t i=0;i<s;++i){
					NT sum=0;
					for(size_t j=0;j<s;++j)
						sum+=Ai[i][j]*dif[j];
					A1.push_back(sum);
        }
        // 				
        // A1 = A1/1+e_i^T*A^{-1}*dif
        NT Ipp1(1+A1[pos]);
        for(size_t i=0;i<s;++i)
					A1[i]/=Ipp1;
				// 
				// Api means "A prime inverse"
				// Api = A^{-1} - A1 * A^{-1}[pos]
				matrix_type Api; 
        tvector<NT> r;
        Api.reserve(s);
        r.reserve(s);
        for(size_t i=0;i<s;++i){
					r.clear();
					for(size_t j=0;j<s;++j)
						r.push_back(Ai[i][j]-A1[i]*Ai[pos][j]);
					Api.push_back(r);
        }
        return DetData(Api,Ipp1*dA,colsAp);
};
 
} // namespace HeaDaCHe

#endif // DYNAMIC_DETERMINANT_H
// vim: ts=2:expandtab
