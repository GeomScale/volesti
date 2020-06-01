#ifndef MISC_LP_H
#define MISC_LP_H


#include <stdio.h>
#include <cmath>
#include <exception>
#undef Realloc
#undef Free
#include "lp_lib.h"


// compute the chebychev ball of an H-polytope described by a dxd matrix A and  d-dimensional vector b, s.t.: Ax<=b
template <typename NT, class Point, class MT, class VT>
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
    free(row);
    free(colno);

    return res;
}



#endif
