// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#ifndef POINT_H
#define POINT_H


template <class K>
class point
{
private:
    int d;
    typedef std::vector<typename K::FT> Coeff;
    Coeff coeffs;
    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    typedef typename K::FT 	FT;

    point() {}
    
    point(const int dim){
        d=dim;
        coeffs=Coeff(d,0);
    }
    
    point(const int dim, iter begin, iter end){
        d=dim;
        coeffs=Coeff(begin,end);
    }
    
    int dimension(){
        return d;
    }
    
    void set_dimension(const int dim){
        d=dim;
    }
    
    void set_coord(const int i, FT coord){
        coeffs[i]=coord;
    }
    
    FT operator[] (const int i){
        return coeffs[i];
    }
    
    point operator+ (point& p){
        point temp(p.dimension());
        typename Coeff::iterator tmit=temp.iter_begin();
        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit, ++tmit){
			(*tmit)=(*pit)+(*mit);
		}
        //for (int i=0; i<d; i++){
            //temp.coeffs[i]=p[i]+coeffs[i];
        //}
        return temp;
    }
    
    point operator- (point& p){
        point temp(p.dimension());
        typename Coeff::iterator tmit=temp.iter_begin();
        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit, ++tmit){
			(*tmit)=(*mit)-(*pit);
		}
        //for (int i=0; i<d; i++){
            //temp.coeffs[i]=coeffs[i]-p[i];
        //}
        
        return temp;
    }
    
    point operator* (const FT& k){
        point temp(d,iter_begin(),iter_end());
        //point temp=*this;
        typename Coeff::iterator tmit=temp.iter_begin();
        for( ; tmit<temp.iter_end(); ++tmit){
			(*tmit)*=k;
		}
        //for (int i=0; i<d; i++){
          //  temp.coeffs[i]=coeffs[i]*k;
        //}
        
        return temp;
    }


    FT dot(point& p){
        FT res=FT(0);

        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit){
            res+=(*mit)*(*pit);
        }
        return res;
    }
    
    
    FT squared_length(){
        
        FT lsq=FT(0);
        
        for (int i=0; i<d; i++){
            lsq+=coeffs[i]*coeffs[i];
        }    
        return lsq;
    }
    
    
    iter iter_begin(){
        
        return coeffs.begin();
        
    }
    
    iter iter_end(){
        
        return coeffs.end();
        
    }
    
    
};

template<class K>
point<K> operator* (const typename K::FT& k, point<K>& p){
        
        return p*k;
        
}

#endif
