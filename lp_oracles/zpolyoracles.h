// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

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


#ifndef ZPOLYORACLES_H
#define ZPOLYORACLES_H


#include <stdio.h>
#include <cmath>
#include <exception>
#undef Realloc
#undef Free
#include "lp_lib.h"


template <typename MT, typename Point, typename NT>
bool memLP_Zonotope(const MT &V, const Point &q, NT *row, int *colno){

    //typedef typename Point::FT NT;
    int d=q.dimension();
    lprec *lp;
    int Ncol=V.rows(), j, i;

    try
    {
        lp = make_lp(d, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    //try
    //{
    //colno = (int *) malloc(Ncol * sizeof(*colno));
    //row = (REAL *) malloc(Ncol * sizeof(*row));
    //}
    //catch (std::exception &e)
    //{
    //#ifdef VOLESTI_DEBUG
    //std::cout<<"Linear Program for membership failed "<<e.what()<<std::endl;
    //#endif
    //return false;
    //}

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i = 0;  i< d; ++i) {
        /* construct all rows */
        for(j=0; j<Ncol; j++){
            colno[j] = j+1;
            row[j] = V(j,i);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, q[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the bounds
    for(j=0; j<Ncol; j++){
        //colno[j] = j+1; /* j_th column */
        row[j] = 0.0;
        set_bounds(lp, j+1, -1.0, 1.0);
    }

    // set the objective function
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
        return false;
    }
    delete_lp(lp);
    return true;
}


// compute the intersection of a ray with a V-polytope
// if maxi is true compute positive lambda, when the ray is p + lambda \cdot v
// otherwise compute the negative lambda
template <typename NT, typename MT, typename Point>
std::pair<NT,NT> intersect_line_zono(const MT &V, const Point &p, const Point &v, NT *row, int *colno){

    std::pair<NT,NT> pair_res;
    int d=v.dimension(), i;
    lprec *lp;//, *lp2;
    int m=V.rows();
    m++;
    int Ncol=m, j, Nrows;
    NT res;
    Nrows = d;

    try
    {
        lp = make_lp(Nrows, Ncol);
        //lp2 = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        // return -1.0;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
    // set_add_rowmode(lp2, TRUE);
    for (i=0; i<d; i++){
        /* construct all rows  */
        for(j=0; j<m-1; j++){
            colno[j] = j+1; /* j_th column */
            row[j] = V(j,i);
        }
        colno[m-1] = m; /* last column */
        row[m-1] = v[i];

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, m, row, colno, EQ, p[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        }

    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<m-1; j++){
        //colno[j] = j+1; /* j_th column */
        set_bounds(lp, j + 1, -1.0, 1.0);
        row[j] = 0;
    }
    //colno[m - 1] =m; /* last column */
    row[m-1] = 1.0;
    set_bounds(lp, m, -infinite, infinite);

    // set objective function
    try
    {
        if(!set_obj_fnex(lp, m, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
    }

    //int* bas = (int *)malloc((d+m+1) * sizeof(int));
    set_maxim(lp);
    set_verbose(lp, NEUTRAL);
    solve(lp);
    pair_res.second = NT(-get_objective(lp));
    set_minim(lp);
    solve(lp);
    pair_res.first = NT(-get_objective(lp));

    delete_lp(lp);
    return pair_res;
}


#endif
