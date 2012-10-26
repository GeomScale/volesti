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

#ifndef EIGEN_FUNCTIONS_H
#define EIGEN_FUNCTIONS_H

#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <Eigen/Eigen>

namespace HeaDaCHe {

typedef CGAL::Linear_algebraCd<CGAL::Gmpz>      LA;
typedef LA::Matrix                     					LAMatrixZ;

template <class _Matrix>
typename _Matrix::NT eigen_determinant(const _Matrix &m){
        typedef _Matrix                                 Matrix;
        typedef typename Matrix::NT                     NT;
				typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>
				                                                EigMatrix;
        int n=m.column_dimension();
        EigMatrix mat(n,n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        mat(i,j)=m(i,j);
        Eigen::FullPivLU<EigMatrix> lu(mat);
        return lu.determinant();
        //return mat.determinant();
}

template < >
CGAL::Gmpz eigen_determinant(const LAMatrixZ &m){
        typedef Eigen::Matrix<CGAL::Gmpq,Eigen::Dynamic,Eigen::Dynamic>
				                                                EigMatrix;
        int n=m.column_dimension();
        EigMatrix mat(n,n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        mat(i,j)=m(i,j);
        Eigen::FullPivLU<EigMatrix> lu(mat);
        CGAL_assertion(lu.determinant().denominator() == 1);
        return CGAL::Gmpz(lu.determinant().numerator());
        //return 
        //return mat.determinant();
}

template <class _Matrix>
_Matrix eigen_inverse(const _Matrix &m){
        typedef _Matrix                                 Matrix;
        typedef typename Matrix::NT                     NT;
        typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>
				                                                EigMatrix;

        int n=m.column_dimension();
        EigMatrix mat(n,n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        mat(i,j)=m(i,j);
        mat=mat.inverse();
        Matrix inv(n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        inv(i,j)=mat(i,j);
        return inv;
}

// bug: incorrect results
template <class _Matrix>
_Matrix eigen_adjoint(const _Matrix &m){
        typedef _Matrix                                 Matrix;
        typedef typename Matrix::NT                     NT;
        typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>
				                                                EigMatrix;
        int n=m.column_dimension();
        EigMatrix mat(n,n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        mat(i,j)=m(i,j);
        mat.adjointInPlace();
        //mat=mat.adjoint();
        Matrix adj(n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        adj(i,j)=mat(i,j);
        return adj;
}

template <class _Matrix>
_Matrix eigen_naive_adjoint(const _Matrix &m){
        typedef _Matrix                                 Matrix;
        typedef typename Matrix::NT                     NT;
        typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>
				                                                EigMatrix;
				                                             
        int n=m.column_dimension();
        EigMatrix mat(n,n);
        EigMatrix adjminor(n-1,n-1);
        /*
        EigMatrix m1(3,3);
        m1(0,0)=1; m1(0,1)=2; m1(0,2)=3;
        m1(1,0)=5; m1(1,1)=5; m1(1,2)=6;
        m1(2,0)=7; m1(2,1)=8; m1(2,2)=9;
        
        //Eigen::FullPivLU<EigMatrix> lu(m1);
        std::cout<<"simple_det="<<m1.determinant()<<std::endl;
        */
        for(size_t i=0;i<n;++i)
					for(size_t j=0;j<n;++j){
						for(size_t im=0;im<n;++im)
						  for(size_t jm=0;jm<n;++jm){
								//std::cout << i << " "<< j << " "<< im << " "<< jm<<std::endl;
								if(im<i && jm<j)
								  adjminor(im,jm)=m(im,jm);
								else if(im<i && jm>j)
									adjminor(im,jm-1)=m(im,jm);
								else if(im>i && jm<j)
									adjminor(im-1,jm)=m(im,jm);
								else if(im>i && jm>j)
									adjminor(im-1,jm-1)=m(im,jm);
								//else std::cout << "no"<<std::endl;
							}
						NT sgn = (i+j)%2 ? -1 : 1;
						Eigen::FullPivLU<EigMatrix> lu(adjminor);
						mat(j,i)=sgn*lu.determinant();					
						}
        Matrix adj(n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j)
                        adj(i,j)=mat(i,j);
        return adj;
}

template < >
LAMatrixZ eigen_naive_adjoint(const LAMatrixZ &m){
        typedef Eigen::Matrix<CGAL::Gmpq,Eigen::Dynamic,Eigen::Dynamic>
				                                                EigMatrix;
				
        int n=m.column_dimension();
        EigMatrix mat(n,n);
        EigMatrix adjminor(n-1,n-1);
        
        EigMatrix m1(3,3);
        m1(0,0)=1; m1(0,1)=2; m1(0,2)=7;
        m1(1,0)=5; m1(1,1)=5; m1(1,2)=6;
        m1(2,0)=7; m1(2,1)=8; m1(2,2)=9;
        
        //Eigen::FullPivLU<EigMatrix> lu(m1);
        //std::cout<<"SP_simple_det="<<m1.determinant()<<std::endl;
        //exit(1);
        for(size_t i=0;i<n;++i)
					for(size_t j=0;j<n;++j){
						for(size_t im=0;im<n;++im)
						  for(size_t jm=0;jm<n;++jm){
								//std::cout << i << " "<< j << " "<< im << " "<< jm<<std::endl;
								if(im<i && jm<j)
								  adjminor(im,jm)=m(im,jm);
								else if(im<i && jm>j)
									adjminor(im,jm-1)=m(im,jm);
								else if(im>i && jm<j)
									adjminor(im-1,jm)=m(im,jm);
								else if(im>i && jm>j)
									adjminor(im-1,jm-1)=m(im,jm);
								//else std::cout << "no"<<std::endl;
							}
						CGAL::Gmpq sgn = (i+j)%2 ? -1 : 1;
						Eigen::FullPivLU<EigMatrix> lu(adjminor);
						mat(j,i)=sgn*lu.determinant();					
						}
        LAMatrixZ adj(n);
        for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j){
                        adj(i,j)=mat(i,j).numerator();
                        CGAL_assertion(mat(i,j).denominator() == 1);
                }
        return adj;
}

} // namespace HeaDaCHe

#endif // EIGEN_FUNCTIONS_H
