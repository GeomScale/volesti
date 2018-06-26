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
    

//public:
    
    typedef typename std::vector<K>::iterator iter;
    
    point() {}
    
    point(int dim){
        d=dim;
        coeffs=coeff(d,0);
        //for (int i=0; i<d; i++){
            //coeffs.push_back(0.0);
           // coeffs[i]=0.0;
        //}
    }
    
    point(int dim, iter begin, iter end){
        d=dim;
        coeffs=coeff(begin,end);
        //for (iter it=begin; it!=end; ++it){
            //coeffs.push_back(*it);
        //}
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
        //temp.set_dimension(d);
        for (int i=0; i<d; i++){
            temp.coeffs[i]=p[i]+coeffs[i];
        }
        
        return temp;
    }
    
    point operator- (point& p){
        point temp(d);
        //temp.set_dimension(d);
        for (int i=0; i<d; i++){
            //temp.coeffs.push_back(coeffs[i]-p[i]);
            temp.coeffs[i]=coeffs[i]-p[i];
        }
        
        return temp;
    }
    
    point operator* (const K& k){
        point temp(d);
        //temp.set_dimension(d);
        for (int i=0; i<d; i++){
            //temp.coeffs.push_back(coeffs[i]*k);
            temp.coeffs[i]=coeffs[i]*k;
        }
        
        return temp;
    }
    
    
    K squared_length(){
        
        K lsq=K(0);
        
        for (int i=0; i<d; i++){
            lsq+=coeffs[i]*coeffs[i];
        }
        //lsq=std::sqrt(lsq);
        
    
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
