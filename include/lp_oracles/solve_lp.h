// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

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


#ifndef SOLVE_LP_H
#define SOLVE_LP_H


#include <stdio.h>
#include <cmath>
#include <exception>
#undef Realloc
#undef Free
#include "lp_lib.h"


// compute the chebychev ball of an H-polytope described by a dxd matrix A and  d-dimensional vector b, s.t.: Ax<=b
template <typename NT, typename Point, typename MT, typename VT>
std::pair<Point,NT> ComputeChebychevBall(MT &A, VT &b){

    lprec *lp;
    int d = A.cols();
    int Ncol=d+1, j, m=A.rows(), i;
    int *colno = NULL;

    REAL *row = NULL;
    std::pair<Point,NT> exception_pair(Point(1),-1.0);

    try
    {
        lp = make_lp(m, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for chebychev center "<<e<<std::endl;
        #endif
        return exception_pair;
    }
    
    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    /* create space large enough for one row */
    try
    {
        colno = (int *) malloc(Ncol * sizeof(*colno));
        row = (REAL *) malloc(Ncol * sizeof(*row));
    }
    catch (std::exception &e)
    {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Linear Program for chebychev center failed "<<e.what()<<std::endl;
        #endif
        return exception_pair;
    }

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    NT sum;
    for (i = 0; i < m; ++i) {
        /* construct all rows */
        sum=NT(0);
        for(j=0; j<d; j++){
            colno[j] = j+1;
            row[j] = A(i,j);
            sum+=A(i,j)*A(i,j);
        }
        colno[d] = d+1; /* last column */
        row[d] = std::sqrt(sum);

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, d+1, row, colno, LE, b(i))) throw false;
        }
        catch (bool e)
        {
            #ifdef VOLESTI_DEBUG
            std::cout<<"Could not define constriants for the Linear Program for chebychev center "<<e<<std::endl;
            #endif
            return exception_pair;
        }
    }

    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
    for (j = 0; j < d; j++) {
        colno[j] = j + 1;
        row[j] = 0;
        set_bounds(lp, j + 1, -infinite, infinite);
    }
    colno[d] = d + 1; /* last column */
    row[d] = 1.0;
    set_bounds(lp, d + 1, 0.0, infinite);

	// set the objective function
    try
    {
        if (!set_obj_fnex(lp, d + 1, row, colno)) throw false;
    }
    catch (bool e)
    {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Could not define objective function for the Linear Program for chebychev center "<<e<<std::endl;
        #endif
        return exception_pair;
    }

    /* set the object direction to maximize */
    set_maxim(lp);

    /* I only want to see important messages on screen while solving */
    set_verbose(lp, NEUTRAL);

    /* Now let lpsolve calculate a solution */
    try
    {
        if (solve(lp) != OPTIMAL) throw false;
    }
    catch (bool e)
    {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Could not solve the Linear Program for chebychev center "<<e<<std::endl;
        #endif
        return exception_pair;
    }

    std::pair<Point,NT> res;

    std::vector<NT> temp_p(d,0);
    get_variables(lp, row);
    for(j = 0; j < d; j++){
        temp_p[j]=NT(row[j]);
    }
    Point xc( d , temp_p.begin() , temp_p.end() );
    NT r=NT(get_objective(lp));
    res = std::pair<Point,NT> (xc,r);
    delete_lp(lp);

    return res;
}


template <typename VT, typename MT, typename Point>
Point PointInIntersection(MT V1, MT V2, Point direction, bool &empty) {

    typedef typename Point::FT NT;
    unsigned int d = V1.cols();
    unsigned int k1 = V1.rows();
    unsigned int k2 = V2.rows();
    unsigned int k = k1 + k2;
    VT cb(k1);
    lprec *lp;
    int Ncol=k, *colno = NULL, j, i;
    REAL *row = NULL;
    Point p(d);

    try
    {
        lp = make_lp(d+2, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    try
    {
        colno = (int *) malloc(Ncol * sizeof(*colno));
        row = (REAL *) malloc(Ncol * sizeof(*row));
    }
    catch (std::exception &e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Linear Program for membership failed "<<e.what()<<std::endl;
#endif
        return false;
    }

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */


    for (i = 0;  i< d+2; ++i) {
        /* construct all rows */
        for(j=0; j<k1; j++){
            colno[j] = j+1;
            if (i==d) {
                row[j] = 1.0;
            } else if(i==d+1){
                row[j] = 0.0;
            } else {
                row[j] = V1(j, i);
            }
        }
        for(j=0; j<k2; j++){
            colno[k1+j] = k1+j+1;
            if (i==d) {
                row[k1+j] = 0.0;
            } else if(i==d+1){
                row[k1+j] = 1.0;
            } else {
                row[k1+j] = -V2(j, i);
            }
        }

        /* add the row to lpsolve */
        try {
            if(i==d || i==d+1) {
                if (!add_constraintex(lp, Ncol, row, colno, EQ, 1.0)) throw false;
            } else {
                if (!add_constraintex(lp, Ncol, row, colno, EQ, 0.0)) throw false;
            }
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }

    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    const NT* direction_data = direction.getCoefficients().data();
    REAL* row_temp = row;

    // set the bounds
    for(j=0; j<Ncol; ++j){
        colno[j] = j+1; /* j_th column */
        *row_temp = *direction_data;
        row_temp++;
        direction_data++;
        set_bounds(lp, j+1, 0.0, infinite);
    }

    try
    {
        if(!set_obj_fnex(lp, Ncol, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    /* set the object direction to maximize */
    set_maxim(lp);

    /* I only want to see important messages on screen while solving */
    set_verbose(lp, NEUTRAL);

    /* Now let lpsolve calculate a solution */
    if (solve(lp) != OPTIMAL){
        delete_lp(lp);
        empty = true;
        return p;
    }
    get_variables(lp, row);
    delete_lp(lp);

    for ( j=0; j<k1; ++j) {
        cb(j) = row[j];
    }

    p = V1.transpose()*cb;
    empty = false;
    return p;

}


#endif
