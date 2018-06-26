#ifndef POINT_H
#define POINT_H


template <class Cr>
class point
{
public:
    typedef typename Cr::RT 	K;
    int d;
    typedef std::vector<K> coeff;
    coeff coeffs;
    
    typedef typename std::vector<K>::iterator iter;
    
    point() {}
    
    point(int dim){
        d=dim;
        coeffs=coeff(d,0);
    }
    
    point(int dim, iter begin, iter end){
        d=dim;
        coeffs=coeff(begin,end);
    }
    
    int dimension(){
        return d;
    }
    
    void set_dimension(int dim){
        d=dim;
    }
    
    void set_coord(int i, K coord){
        coeffs[i]=coord;
    }
    
    K operator[] (const int i){
        return coeffs[i];
    }
    
    point operator+ (point& p){
        point temp(d);
        for (int i=0; i<d; i++){
            temp.coeffs[i]=p[i]+coeffs[i];
        }
        
        return temp;
    }
    
    point operator- (point& p){
        point temp(d);
        for (int i=0; i<d; i++){
            temp.coeffs[i]=coeffs[i]-p[i];
        }
        
        return temp;
    }
    
    point operator* (const K& k){
        point temp(d);
        for (int i=0; i<d; i++){
            temp.coeffs[i]=coeffs[i]*k;
        }
        
        return temp;
    }
    
    
    K squared_length(){
        
        K lsq=K(0);
        
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

template<class Cr>
point<Cr> operator* (const typename Cr::RT& k, point<Cr>& p){
        
        return p*k;
        
}

#endif
