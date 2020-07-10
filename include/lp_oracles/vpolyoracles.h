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


#ifndef VPOLYORACLES_H
#define VPOLYORACLES_H


#include <stdio.h>
#include <cmath>
#include <exception>
#undef Realloc
#undef Free
#include "lp_lib.h"


// return true if q belongs to the convex hull of the V-polytope described by matrix V
// otherwise return false
template <typename MT, typename Point, typename NT>
bool memLP_Vpoly(const MT &V, const Point &q, NT *row, int *colno){

    //typedef typename Point::FT NT;
    int d=q.dimension();
    lprec *lp;
    int Ncol=d+1, j, i, m=V.rows();
    m++;

    try
    {
        lp = make_lp(m, Ncol);
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
    // std::cout<<"Linear Program for membership failed "<<e.what()<<std::endl;
    //#endif
    //return false;
    //}

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i = 0;  i< m-1; ++i) {
        /* construct all rows */
        for(j=0; j<d; j++){
            colno[j] = j+1;
            row[j] = V(i,j);
        }
        colno[d] = d+1;
        row[d] = -1.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, d+1, row, colno, LE, 0.0)) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }
    for(j=0; j<d; j++){
        //colno[j] = j+1; // last column
        row[j] = q[j];
    }
    //colno[d] = d+1; // last column
    row[d] = -1.0;

    /* add the row to lpsolve */
    try {
        if(!add_constraintex(lp, d+1, row, colno, LE, 1.0)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the bounds
    for(j=0; j<d; j++){
        //colno[j] = j+1; /* j_th column */
        row[j] = q[j];
        set_bounds(lp, j+1, -infinite, infinite);
    }
    //colno[d] = d+1; /* last column */
    row[d] = -1.0;
    set_bounds(lp, d+1, -infinite, infinite);

    // set the objective function
    try
    {
        if(!set_obj_fnex(lp, d+1, row, colno)) throw false;
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
    try
    {
        if (solve(lp) != OPTIMAL) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not solve the Linear Program for memebrship "<<e<<std::endl;
#endif
        return false;
    }

    NT r = NT(get_objective(lp));
    delete_lp(lp);
    if(r>0.0){
        return false;
    }
    return true;
}



// compute the intersection of a ray with a V-polytope
// if maxi is true compute positive lambda, when the ray is p + lambda \cdot v
// otherwise compute the negative lambda
template <typename NT, typename MT, typename Point>
NT intersect_line_Vpoly(const MT &V, const Point &p, const Point &v,
        NT *conv_comb, NT *row, int *colno,  bool maxi, bool zonotope){

    int d=v.dimension(), i;
    lprec *lp;
    int m=V.rows();
    m++;
    int Ncol=m, j, Nrows;
    NT res;
    if(!zonotope) {
        Nrows = d+1;
    } else {
        Nrows = d;
    }

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return -1.0;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

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
            return -1.0;
        }

    }

    if(!zonotope) {
        for (j = 0; j < m - 1; j++) {
            //colno[j] = j + 1; /* j_th column */
            row[j] = 1.0;
        }
        //colno[m - 1] = m; /* last column */
        row[m - 1] = 0.0;

        /* add the row to lpsolve */
        try {
            if (!add_constraintex(lp, m, row, colno, EQ, 1.0)) throw false;
        }
        catch (bool e) {
#ifdef VOLESTI_DEBUG
            std::cout << "Could not construct constaints for the Linear Program for ray-shooting " << e << std::endl;
#endif
            return -1.0;
        }
    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<m-1; j++){
        //colno[j] = j+1; /* j_th column */
        if (!zonotope) {
            set_bounds(lp, j + 1, 0.0, 1.0);
        } else {
            set_bounds(lp, j + 1, -1.0, 1.0);
        }
        row[j] = 0;
    }
    //colno[m - 1] = m; /* last column */
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
        return -1.0;
    }

    if(maxi) {  /* set the object direction to maximize */
        set_maxim(lp);
    }else{      /* set the object direction to minimize */
        set_minim(lp);
    }
    set_verbose(lp, NEUTRAL);

    /* Now let lpsolve calculate a solution */
    try
    {
        if (solve(lp) != OPTIMAL) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not solve the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return -1.0;
    }

    res = NT(-get_objective(lp));
    //std::cout<<"res = "<<res<<std::endl;
    get_variables(lp, conv_comb);
    delete_lp(lp);
    return res;
}


template <typename NT, typename MT, typename Point>
std::pair<NT,NT> intersect_double_line_Vpoly(const MT &V, const Point &p, const Point &v, NT *row, int *colno){

    int d=v.dimension(), i;
    lprec *lp;
    int m=V.rows();
    m++;
    int Ncol=m, j, Nrows;
    NT res;
    Nrows = d+1;
    std::pair<NT,NT> res_pair;

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return res_pair;
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
    //    std::cout<<"Linear Program for ray-shooting failed "<<e.what()<<std::endl;
//#endif
    //      return res_pair;
    //}

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

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
            return res_pair;
        }

    }

    for (j = 0; j < m - 1; j++) {
        colno[j] = j + 1; /* j_th column */
        row[j] = 1.0;
    }
    colno[m - 1] = m; /* last column */
    row[m - 1] = 0.0;

    /* add the row to lpsolve */
    try {
        if (!add_constraintex(lp, m, row, colno, EQ, 1.0)) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout << "Could not construct constaints for the Linear Program for ray-shooting " << e << std::endl;
#endif
        return res_pair;
    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<m-1; j++){
        colno[j] = j+1; /* j_th column */
        set_bounds(lp, j + 1, 0.0, 1.0);
        row[j] = 0;
    }
    colno[m - 1] =m; /* last column */
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
        return res_pair;
    }

    set_maxim(lp);
    set_verbose(lp, NEUTRAL);
    solve(lp);

    res_pair.second = NT(-get_objective(lp));

    set_minim(lp);
    solve(lp);
    res_pair.first = NT(-get_objective(lp));


    delete_lp(lp);
    return res_pair;
}


#endif
