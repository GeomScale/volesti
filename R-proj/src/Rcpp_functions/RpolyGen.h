// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RPOLYGEN_H
#define RPOLYGEN_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix poly_gen (int kind_gen, bool Vpoly_gen, int dim_gen, int m_gen){

    Rcpp::NumericMatrix Mat;
    if (kind_gen == 0) {
        Zonotope ZP = gen_zonotope<Zonotope, RNGType>(dim_gen, m_gen);
        Mat = extractMatPoly(ZP);
    } else if (Vpoly_gen) {
        if (kind_gen == 1) {
            VP = gen_cube<Vpolytope>(dim_gen, true);
            Mat = extractMatPoly(VP);
        } else if (kind_gen == 2) {
            VP = gen_cross<Vpolytope>(dim_gen, true);
            Mat = extractMatPoly(VP);
        } else if (kind_gen == 3) {
            VP = gen_simplex<Vpolytope>(dim_gen, true);
            Mat = extractMatPoly(VP);
        } else {
            return Mat;
        }
    } else {
        if (kind_gen == 1) {
            HP = gen_cube<Hpolytope>(dim_gen, false);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 2) {
            HP = gen_cross<Hpolytope>(dim_gen, false);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 3) {
            HP = gen_simplex<Hpolytope>(dim_gen, false);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 4) {
            HP = gen_prod_simplex<Hpolytope>(dim_gen);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 5) {
            HP = gen_skinny_cube<Hpolytope>(dim_gen);
            Mat = extractMatPoly(HP);
        } else {
            return Mat;
        }
    }
    return Mat;

}



#endif
