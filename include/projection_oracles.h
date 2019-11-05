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


#ifndef PROJECTION_ORACLES_H
#define PROJECTION_ORACLES_H

template <class MT, class VT, class Point, typename NT>
NT intersect_line_proj_poly(MT &T, MT &A, VT &b, Point &r, Point &v,  NT *conv_comb, NT *row, int *colno) {

    int d = v.dimension(), m = A.rows(), i, Ncol = T.cols()+1, j, Nrows, k = T.cols();
    lprec *lp;
    NT res;
    Nrows = m + d;

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
        for(j=0; j<k; j++){
            colno[j] = j+1; /* j_th column */
            row[j] = T(i,j);
        }
        colno[k] = k + 1; /* last column */
        row[k] = v[i];

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, r[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return -1.0;
        }

    }

    for (i=0; i<m; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            //colno[j] = j+1; /* j_th column */
            row[j] = A(i,j);
        }
        //colno[k] = k + 1; /* last column */
        row[k] = 0.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, LE, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return -1.0;
        }

    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<Ncol; j++){
        set_bounds(lp, j + 1, -infinite, infinite);
        //colno[j] = j+1; /* j_th column */
        row[j] = 0;
    }
    row[Ncol - 1] = 1.0;

    // set objective function
    try
    {
        if(!set_obj_fnex(lp, Ncol, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return -1.0;
    }

        /* set the object direction to minimize*/
    set_minim(lp);
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
    std::cout<<"res = "<<res<<std::endl;
    get_variables(lp, conv_comb);
    delete_lp(lp);
    return res;

}


template <class MT, class VT, class Point, typename NT>
std::pair<NT,NT> intersect_double_line_proj_poly(MT &T, MT &A, VT &b, Point &r, Point &v,  NT *conv_comb, NT *row, int *colno) {

    int d = v.dimension(), m = A.rows(), i, Ncol = T.cols()+1, j, Nrows, k = T.cols();
    lprec *lp;
    NT res;
    Nrows = m + d;
    std::pair<NT,NT> exc_res(0,0);

    try
    {
        lp = make_lp(Nrows, Ncol);
        if(lp == NULL) throw false;
    }
    catch (bool e) {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return exc_res;
    }

    REAL infinite = get_infinite(lp); /* will return 1.0e30 */

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    for (i=0; i<d; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            colno[j] = j+1; /* j_th column */
            row[j] = T(i,j);
        }
        colno[k] = k + 1; /* last column */
        row[k] = v[i];

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, EQ, r[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    for (i=0; i<m; i++){
        /* construct all rows  */
        for(j=0; j<k; j++){
            //colno[j] = j+1; /* j_th column */
            row[j] = A(i,j);
        }
        //colno[k] = k + 1; /* last column */
        row[k] = 0.0;

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, Ncol, row, colno, LE, b(i))) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
            return exc_res;
        }

    }

    //set the bounds
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // set the objective function
    for(j=0; j<Ncol; j++){
        set_bounds(lp, j + 1, -infinite, infinite);
        //colno[j] = j+1; /* j_th column */
        row[j] = 0;
    }
    row[Ncol - 1] = 1.0;

    // set objective function
    try
    {
        if(!set_obj_fnex(lp, Ncol, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for ray-shooting "<<e<<std::endl;
#endif
        return exc_res;
    }

    std::pair<NT,NT> res_pair;
    set_verbose(lp, NEUTRAL);

    set_maxim(lp);
    solve(lp);

    res_pair.second = NT(-get_objective(lp));

    set_minim(lp);
    solve(lp);
    res_pair.first = NT(-get_objective(lp));

    delete_lp(lp);

    //std::cout<<"l1 = "<<res_pair.first<<" l2 = "<<res_pair.second<<std::endl;
    return res_pair;

}



template <class MT, class VT, class Point>
bool memLP_proj_poly(MT &T, MT &A, VT &b, Point &q) {

    typedef typename Point::FT NT;
    int d=q.dimension(), m = A.rows(), k = T.cols();
    lprec *lp;
    int Ncol=k, *colno = NULL, j, i;

    REAL *row = NULL;

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

    for (i = 0;  i< d; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            colno[j] = j+1;
            row[j] = T(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, row, colno, EQ, q[i])) throw false;
        }
        catch (bool e)
        {
#ifdef VOLESTI_DEBUG
            std::cout<<"Could not construct constaints for the Linear Program for membership "<<e<<std::endl;
#endif
            return false;
        }
    }

    for (i = 0;  i< m; ++i) {
        /* construct all rows */
        for(j=0; j<k; j++){
            //colno[j] = j+1;
            row[j] = A(i,j);
        }

        /* add the row to lpsolve */
        try {
            if(!add_constraintex(lp, k, row, colno, LE, b(i))) throw false;
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

    for(j=0; j<k; j++){
        //colno[j] = j+1; // last column
        row[j] = 0.0;
        set_bounds(lp, j+1, -infinite, infinite);
    }

    /* add the row to lpsolve */
    try {
        if(!set_obj_fnex(lp, k, row, colno)) throw false;
    }
    catch (bool e)
    {
#ifdef VOLESTI_DEBUG
        std::cout<<"Could not construct objective function for the Linear Program for membership "<<e<<std::endl;
#endif
        return false;
    }

    //set the bounds

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


#endif
