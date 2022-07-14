
#include <fstream>

template <typename Point>
void write_to_file(std::string filename, std::vector<Point> const& randPoints) {
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); //save and redirect
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf); //reset to standard output again
}

template<typename MT>
MT rand_matrix(int n){
    int i,j;
    double val;
    MT A = MT::Identity(n,n);
    for(i = 0; i < n; i++){
        for(j = i+1; j < n; j++){
            val = 2* ((double) rand()/RAND_MAX) - 1;
            A(j,i) = val;
        }
    }
    return A;
}

template<typename Point, typename RNGType>
Point getDirection(unsigned int const& dim, RNGType &rng, bool normalize=true){
    double normal = 0.;
    Point p(dim);
    double* data = p.pointerToData();

    for (unsigned int i=0; i<dim; ++i){
        *data = rng.sample_ndist();
        normal += *data * *data;
        data++;
    }

    normal = 1./std::sqrt(normal);
    if (normalize) p *= normal;
    return p;
}

template<typename VT, typename MT>
MT rebuildMatrix(const VT &xvector, const unsigned int n){
    MT A = MT::Identity(n,n);
    double coeff;
    for(int i = 0; i < n ; ++i){
        for(int j = i+1; j < n; ++j){
            int ind = ((((n<<1)-i-2)*(i+1)) >> 1)  + j - n;
            coeff = xvector[ind];
            A(i,j) = coeff;
            A(j,i) = coeff;
        }
    }
    return A;
}

template<typename VT, typename MT>
bool membership(const VT &xvector, const unsigned int n){
    int d = n*(n-1)/2;
    for(int i = 0; i < d; ++i)
        if((xvector(i) > 1) || (xvector(i) < -1)) return false;
    MT A = rebuildMatrix<VT, MT>(xvector, n);
    Eigen::LDLT<MT> A_ldlt(A);
    if (A_ldlt.info() != Eigen::NumericalIssue && A_ldlt.isPositive())
        return true;
    return false;
}

template<typename VT, typename MT, typename PointList>
void check_output(PointList &randPoints, int num_points, int n){
    for(int i = 0; i < num_points ; ++i){
        if(!membership<VT,MT>(randPoints[i].getCoefficients(), n)){
            std::cout << "ALERT\n";
        }else{ std::cout << "OK\n";}
    }
}