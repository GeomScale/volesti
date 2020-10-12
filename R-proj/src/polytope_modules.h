
class Hpolytope {
public:
    Hpolytope() {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b) : A(_A), b(_b), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_A.ncol()), type(1) {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b, double volume) : A(_A), b(_b), vol(volume), dimension(_A.ncol()), type(1) {}
    Rcpp::NumericMatrix A;
    Rcpp::NumericVector b;
    double vol;
    unsigned int dimension;
    int type;
};

class Vpolytope {
public:
    Vpolytope() {}
    Vpolytope(Rcpp::NumericMatrix _V) : V(_V), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_V.ncol()), type(2) {}
    Vpolytope(Rcpp::NumericMatrix _V, double volume) : V(_V), vol(volume), dimension(_V.ncol()), type(2) {}
    Rcpp::NumericMatrix V;
    double vol;
    unsigned int dimension;
    int type;
};

class Zonotope {
public:
    Zonotope() {}
    Zonotope(Rcpp::NumericMatrix _G) : G(_G), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_G.ncol()), type(3) {}
    Zonotope(Rcpp::NumericMatrix _G, double volume) : G(_G), vol(volume), dimension(_G.ncol()), type(3) {}
    Rcpp::NumericMatrix G;
    double vol;
    unsigned int dimension;
    int type;
};

class VpolytopeIntersection {
public:
    VpolytopeIntersection() {}
    VpolytopeIntersection(Rcpp::NumericMatrix _V1, Rcpp::NumericMatrix _V2) : V1(_V1), V2(_V2), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_V1.ncol()), type(4) {}
    VpolytopeIntersection(Rcpp::NumericMatrix _V1, Rcpp::NumericMatrix _V2, double volume) : V1(_V1), V2(_V2), vol(volume), dimension(_V1.ncol()), type(4) {}
    Rcpp::NumericMatrix V1;
    Rcpp::NumericMatrix V2;
    double vol;
    unsigned int dimension;
    int type;
};