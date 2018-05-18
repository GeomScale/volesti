template <typename K>
class point_d{
private:
    std::vector<K> coeffs;
    int d;
    typedef std::vector<K>::iterator veciter;

public:
    point_d() {}
    
    point_d(int dim){
        d=dim;
        int i;
        for (i=0; i<d; i++){
            coeffs.push_back(0.0);
        }
    }
    
    point_d(int dim, veciter begin, veciter end){
        d=dim;
        for (veciter it=begin; it!=end; ++it){
            coeffs.push_back(*it);
        }
    }
    
    int dimension(){
        return d;
    }
    
    void set_coord(int i, K coord){
        coeffs[i]=coord;
    }
    
    K point_d::operator[] (const int i){
        return coeffs[i];
    }
    
    point_d point_d::operator+ (const point_d& point){
        point_d temp;
        int i;
        for (i=0; i<d; i++){
            temp.set_coord(i, coeffs[i]+point[i]);
        }
        
        return temp;
    }
    
    point_d point_d::operator- (const point_d& point){
        point_d temp;
        int i;
        for (i=0; i<d; i++){
            temp.set_coord(i, coeffs[i]-point[i]);
        }
        
        return temp;
    }
    
    point_d point_d::operator* (const K& k){
        point_d temp;
        for (int i=0; i<d; i++){
            temp.set_coord(i, coeffs[i]*k);
        }
        
        return temp;
    }
    
    
};
